/*
######################################################################
#
# moveset.c
#
# Licensed under the GNU General Public License version 3 or later.
#
# Part of the R/ergmgp package
#
# Implementation of the constrained move sets.  See moveset.h.
#
######################################################################
*/

#include "moveset.h"
#include <R.h>
#include <Rmath.h>

/*Build a move set.  rledyads is the R-side encoding of the free dyad set
  (see to_ergm_Cdouble.rlebdm), or NULL if no dyad-level restriction.*/
MoveSet *MoveSetInit(int family, Network *nwp, double *rledyads){
  MoveSet *ms = (MoveSet *) R_Calloc(1, MoveSet);
  ms->family = family;
  ms->directed = nwp->directed_flag;
  ms->n = nwp->nnodes;
  ms->usedyads = (rledyads != NULL);
  if(ms->usedyads){
    double *p = rledyads;
    ms->dyads = unpack_RLEBDM1D(&p);
  }
  ms->eltail = ms->elhead = NULL;
  ms->nedges = 0;
  ms->elalloc = 0;
  MoveSetRefresh(ms, nwp);
  return ms;
}

void MoveSetDestroy(MoveSet *ms){
  if(ms == NULL) return;
  if(ms->eltail) R_Free(ms->eltail);
  if(ms->elhead) R_Free(ms->elhead);
  R_Free(ms);
}

/*Snapshot the edge list.  The degree families enumerate over pairs of
  edges, and we do not want to be walking a live edge tree while toggling
  it, so we take a copy once per event.*/
void MoveSetRefresh(MoveSet *ms, Network *nwp){
  if(ms->family < MS_MULTITOG) return;   /*Dyad families need no snapshot*/
  Edge ne = EDGECOUNT(nwp);
  if(ne > ms->elalloc){
    Edge na = (ne < 64) ? 64 : ne * 2;
    ms->eltail = (Vertex *) R_Realloc(ms->eltail, na, Vertex);
    ms->elhead = (Vertex *) R_Realloc(ms->elhead, na, Vertex);
    ms->elalloc = na;
  }
  Edge i = 0;
  EXEC_THROUGH_NET_EDGES(t, h, e, {
    if(i < ne){ ms->eltail[i] = t; ms->elhead[i] = h; i++; }
  });
  ms->nedges = i;
}

/*Is dyad (t,h) available for toggling, given any dyad-level restriction?*/
static inline Rboolean ms_dyadok(const MoveSet *ms, Vertex t, Vertex h){
  if(!ms->usedyads) return TRUE;
  return GetRLEBDM1D(t, h, &ms->dyads) ? TRUE : FALSE;
}

/*----------------------------------------------------------------------
  Superset sizes.

  For thinning we need the number of candidates we sample from, and we
  need it to be a quantity we can compute in O(1).  Rather than counting
  legal moves (expensive, and state dependent), we sample from a fixed
  product set and assign rate 0 to the illegal members.  Uniformization
  stays exact; illegal draws are simply self-loops.
  ----------------------------------------------------------------------*/
double MoveSetSupersetSize(const MoveSet *ms, Network *nwp){
  double n = (double) ms->n;
  double ne = (double) EDGECOUNT(nwp);
  switch(ms->family){
    case MS_FREE:
      return n*n;                 /*Ordered (t,h) draws, including rejects*/
    case MS_DYAD:
      return (double) ms->dyads.ndyads;
    case MS_ODEGREES:
    case MS_IDEGREES:
      return ne * n;              /*(edge, node) pairs*/
    case MS_DEGREES:
      return ne * ne * 2.0;       /*(edge, edge, re-pairing) triples*/
    case MS_EDGES:
      return ne * n * n;          /*(edge, ordered dyad) pairs*/
  }
  return 0.0;
}

/*----------------------------------------------------------------------
  Uniform sampling from the superset.
  ----------------------------------------------------------------------*/
Rboolean MoveSetSample(const MoveSet *ms, Network *nwp, Move *mv){
  Vertex t, h, k, a1, a2, b1, b2;
  Edge ne = EDGECOUNT(nwp);

  switch(ms->family){
    case MS_FREE: {
      /*Uniform ordered (t,h) draw.  Every move must correspond to exactly one
        element of the superset, so when undirected we *reject* the lower
        triangle rather than folding it onto the upper one -- folding would
        make each dyad twice as likely as the superset size claims, and would
        inflate the total rate accordingly.*/
      t = 1 + (Vertex)(unif_rand() * ms->n);
      h = 1 + (Vertex)(unif_rand() * ms->n);
      if(t == h) return FALSE;
      if(!ms->directed && t > h) return FALSE;
      mv->ntoggles = 1;
      mv->tails[0] = t; mv->heads[0] = h;
      mv->isedge = IS_OUTEDGE(t, h);
      return TRUE;
    }
    case MS_DYAD: {
      if(ms->dyads.ndyads == 0) return FALSE;
      GetRandRLEBDM1D_ITS(&t, &h, &ms->dyads);
      mv->ntoggles = 1;
      mv->tails[0] = t; mv->heads[0] = h;
      mv->isedge = IS_OUTEDGE(t, h);
      return TRUE;
    }
    case MS_ODEGREES:
    case MS_IDEGREES: {
      /*A uniform (edge, node) pair.  The move rewires the edge to (or
        from) that node, keeping the shared endpoint's degree fixed.*/
      if(ne == 0) return FALSE;
      GetRandEdge(&t, &h, nwp);
      k = 1 + (Vertex)(unif_rand() * ms->n);
      if(k == t || k == h) return FALSE;
      if(ms->family == MS_ODEGREES){
        if(IS_OUTEDGE(t, k)) return FALSE;              /*Target already present*/
        if(!ms_dyadok(ms, t, h) || !ms_dyadok(ms, t, k)) return FALSE;
        mv->ntoggles = 2;
        mv->tails[0] = t; mv->heads[0] = h;             /*Dissolve*/
        mv->tails[1] = t; mv->heads[1] = k;             /*Form*/
      }else{
        if(IS_OUTEDGE(k, h)) return FALSE;
        if(!ms_dyadok(ms, t, h) || !ms_dyadok(ms, k, h)) return FALSE;
        mv->ntoggles = 2;
        mv->tails[0] = t; mv->heads[0] = h;
        mv->tails[1] = k; mv->heads[1] = h;
      }
      mv->isedge = 1;   /*The primary dyad is an edge, by construction*/
      return TRUE;
    }
    case MS_EDGES: {
      /*One uniform edge to dissolve, plus one uniform dyad to form.  The
        superset is (edge, ordered dyad), of size ne*n*n; ne is exactly what
        this constraint holds fixed, so the size does not depend on the state,
        which is what uniformization requires.

        As in MS_FREE the dyad draw is ordered, and when undirected we *reject*
        the lower triangle rather than folding it onto the upper one: every
        legal move must correspond to exactly one superset element, or the
        total rate comes out inflated.*/
      if(ne == 0) return FALSE;
      GetRandEdge(&a1, &a2, nwp);
      t = 1 + (Vertex)(unif_rand() * ms->n);
      h = 1 + (Vertex)(unif_rand() * ms->n);
      if(t == h) return FALSE;
      if(!ms->directed && t > h) return FALSE;
      /*The dyad we form must be a non-edge; this also rules out re-forming the
        edge we just dissolved, which would be a no-op.*/
      if(IS_OUTEDGE(t, h)) return FALSE;
      if(!ms_dyadok(ms, a1, a2) || !ms_dyadok(ms, t, h)) return FALSE;
      mv->ntoggles = 2;
      mv->tails[0] = a1; mv->heads[0] = a2;   /*Dissolve*/
      mv->tails[1] = t;  mv->heads[1] = h;    /*Form*/
      mv->isedge = 1;   /*The primary dyad is an edge, by construction*/
      return TRUE;
    }
    case MS_DEGREES: {
      /*Two uniform edges plus one of the two alternative re-pairings.*/
      if(ne < 2) return FALSE;
      GetRandEdge(&a1, &a2, nwp);
      GetRandEdge(&b1, &b2, nwp);
      if(a1 == b1 && a2 == b2) return FALSE;            /*Same edge*/
      if(a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) return FALSE;  /*Shared node*/
      /*Re-pairing {A,B} is the same move as re-pairing {B,A}, so the ordered
        draw generates every move twice.  Keep one canonical ordering, or the
        total rate comes out doubled.  (The edges are node-disjoint here, so
        a1 != b1 and this test is decisive.)*/
      if(a1 > b1) return FALSE;
      Vertex c1, c2, d1, d2;
      if(unif_rand() < 0.5){ c1 = a1; c2 = b2; d1 = b1; d2 = a2; }
      else                 { c1 = a1; c2 = b1; d1 = a2; d2 = b2; }
      /*Undirected: normalise endpoint order*/
      if(c1 > c2){ Vertex tmp = c1; c1 = c2; c2 = tmp; }
      if(d1 > d2){ Vertex tmp = d1; d1 = d2; d2 = tmp; }
      if(IS_OUTEDGE(c1, c2) || IS_OUTEDGE(d1, d2)) return FALSE;
      if(!ms_dyadok(ms, a1, a2) || !ms_dyadok(ms, b1, b2) ||
         !ms_dyadok(ms, c1, c2) || !ms_dyadok(ms, d1, d2)) return FALSE;
      mv->ntoggles = 4;
      mv->tails[0] = a1; mv->heads[0] = a2;   /*Dissolve*/
      mv->tails[1] = b1; mv->heads[1] = b2;   /*Dissolve*/
      mv->tails[2] = c1; mv->heads[2] = c2;   /*Form*/
      mv->tails[3] = d1; mv->heads[3] = d2;   /*Form*/
      mv->isedge = 1;
      return TRUE;
    }
  }
  return FALSE;
}

/*----------------------------------------------------------------------
  Enumeration.
  ----------------------------------------------------------------------*/
void MoveIterInit(MoveIter *it, const MoveSet *ms, Network *nwp){
  it->ms = ms;
  it->nwp = nwp;
  it->d = 0;
  it->hint = 0;
  it->ei = 0;
  it->ej = 1;   /*MS_DEGREES uses unordered pairs ei<ej*/
  it->k = 0;
  it->variant = 0;
  it->done = 0;
}

Rboolean MoveIterNext(MoveIter *it, Move *mv){
  const MoveSet *ms = it->ms;
  Network *nwp = it->nwp;
  Vertex n = ms->n;

  if(it->done) return FALSE;

  switch(ms->family){
    case MS_FREE:
    case MS_DYAD: {
      /*Walk all dyads (skipping the diagonal, and the lower triangle when
        undirected), honouring any free-dyad restriction.*/
      Dyad ndyad = (Dyad)n * (Dyad)n;
      while(TRUE){
        it->d++;
        if(it->d > ndyad){ it->done = 1; return FALSE; }
        Vertex t, h;
        Dyad2TH(&t, &h, it->d, n);
        if(t == h) continue;
        if(!ms->directed && t > h) continue;
        if(ms->family == MS_DYAD && !ms_dyadok(ms, t, h)) continue;
        mv->ntoggles = 1;
        mv->tails[0] = t; mv->heads[0] = h;
        mv->isedge = IS_OUTEDGE(t, h);
        return TRUE;
      }
    }
    case MS_ODEGREES:
    case MS_IDEGREES: {
      /*For each edge, for each node that could receive the rewire.*/
      while(TRUE){
        if(it->ei >= ms->nedges){ it->done = 1; return FALSE; }
        Vertex t = ms->eltail[it->ei], h = ms->elhead[it->ei];
        it->k++;
        if(it->k > n){ it->k = 0; it->ei++; continue; }
        Vertex k = it->k;
        if(k == t || k == h) continue;
        if(ms->family == MS_ODEGREES){
          if(IS_OUTEDGE(t, k)) continue;
          if(!ms_dyadok(ms, t, h) || !ms_dyadok(ms, t, k)) continue;
          mv->ntoggles = 2;
          mv->tails[0] = t; mv->heads[0] = h;
          mv->tails[1] = t; mv->heads[1] = k;
        }else{
          if(IS_OUTEDGE(k, h)) continue;
          if(!ms_dyadok(ms, t, h) || !ms_dyadok(ms, k, h)) continue;
          mv->ntoggles = 2;
          mv->tails[0] = t; mv->heads[0] = h;
          mv->tails[1] = k; mv->heads[1] = h;
        }
        mv->isedge = 1;
        return TRUE;
      }
    }
    case MS_EDGES: {
      /*Every (edge, non-edge) pair: dissolve the one, form the other.  The
        move set has E*(D-E) members, so this is as expensive as MS_DEGREES;
        prefer thinning where the process allows it.*/
      Dyad ndyad = (Dyad)n * (Dyad)n;
      while(TRUE){
        if(it->ei >= ms->nedges){ it->done = 1; return FALSE; }
        Vertex a1 = ms->eltail[it->ei], a2 = ms->elhead[it->ei];
        it->d++;
        if(it->d > ndyad){ it->d = 0; it->ei++; continue; }
        Vertex t, h;
        Dyad2TH(&t, &h, it->d, n);
        if(t == h) continue;
        if(!ms->directed && t > h) continue;
        if(IS_OUTEDGE(t, h)) continue;   /*Must form a non-edge*/
        if(!ms_dyadok(ms, a1, a2) || !ms_dyadok(ms, t, h)) continue;
        mv->ntoggles = 2;
        mv->tails[0] = a1; mv->heads[0] = a2;   /*Dissolve*/
        mv->tails[1] = t;  mv->heads[1] = h;    /*Form*/
        mv->isedge = 1;
        return TRUE;
      }
    }
    case MS_DEGREES: {
      /*Ordered pairs of distinct, node-disjoint edges, times the two
        alternative re-pairings.  We take ei<ej (unordered pairs) because
        the two re-pairings already cover both orientations.*/
      while(TRUE){
        /*Advance the cursor: variant fastest, then ej, then ei.*/
        it->variant++;
        if(it->variant > 2){ it->variant = 1; it->ej++; }
        if(it->ej >= ms->nedges){ it->ei++; it->ej = it->ei + 1; it->variant = 1; }
        if(it->ei + 1 >= ms->nedges || it->ej >= ms->nedges){ it->done = 1; return FALSE; }
        Vertex a1 = ms->eltail[it->ei], a2 = ms->elhead[it->ei];
        Vertex b1 = ms->eltail[it->ej], b2 = ms->elhead[it->ej];
        if(a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;
        Vertex c1, c2, d1, d2;
        if(it->variant == 1){ c1 = a1; c2 = b2; d1 = b1; d2 = a2; }
        else                { c1 = a1; c2 = b1; d1 = a2; d2 = b2; }
        if(c1 > c2){ Vertex tmp = c1; c1 = c2; c2 = tmp; }
        if(d1 > d2){ Vertex tmp = d1; d1 = d2; d2 = tmp; }
        if(IS_OUTEDGE(c1, c2) || IS_OUTEDGE(d1, d2)) continue;
        if(!ms_dyadok(ms, a1, a2) || !ms_dyadok(ms, b1, b2) ||
           !ms_dyadok(ms, c1, c2) || !ms_dyadok(ms, d1, d2)) continue;
        mv->ntoggles = 4;
        mv->tails[0] = a1; mv->heads[0] = a2;
        mv->tails[1] = b1; mv->heads[1] = b2;
        mv->tails[2] = c1; mv->heads[2] = c2;
        mv->tails[3] = d1; mv->heads[3] = d2;
        mv->isedge = 1;
        return TRUE;
      }
    }
  }
  it->done = 1;
  return FALSE;
}
