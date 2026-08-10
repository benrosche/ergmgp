/*
######################################################################
#
# moveset.h
#
# Licensed under the GNU General Public License version 3 or later.
#
# Part of the R/ergmgp package
#
# Move sets for constrained ERGM generating processes.
#
# An unconstrained EGP moves by toggling one dyad.  Under a constraint
# such as ~degrees no single toggle is legal, so a move must become a
# *set* of simultaneous toggles that preserves the constrained
# statistic.  This file defines those move sets and gives two ways to
# work with them:
#
#   enumeration - walk the entire legal move set (MoveIter).  Exact for
#                 every process, but the move set is O(E^2) for ~degrees,
#                 so this is slow.
#   sampling    - draw a candidate uniformly from a fixed *superset* of
#                 the move set whose size does not depend on the current
#                 state (MoveSetSample).  Illegal candidates are simply
#                 rejected.  This is what makes exact uniformization
#                 ("thinning") possible without ever counting the legal
#                 moves.
#
######################################################################
*/
#ifndef EGP_MOVESET_H
#define EGP_MOVESET_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"

/*Move-set families.  These must agree with the EGP_MS_* constants in
  R/constraints.R.*/
#define MS_FREE      0  /*1 toggle, all dyads*/
#define MS_DYAD      1  /*1 toggle, restricted dyad set*/
#define MS_ODEGREES  2  /*2 toggles sharing a tail*/
#define MS_IDEGREES  3  /*2 toggles sharing a head*/
#define MS_DEGREES   4  /*4 toggles (tetrad), undirected*/

#define MOVE_MAXTOG  4

/*A single move: a set of dyads to toggle simultaneously.*/
typedef struct {
  int ntoggles;
  Vertex tails[MOVE_MAXTOG];
  Vertex heads[MOVE_MAXTOG];
  int isedge;   /*Edge state of the first ("primary") dyad before the move*/
} Move;

/*A move set: the family, plus whatever auxiliary state it needs.*/
typedef struct {
  int family;
  int directed;
  Vertex n;
  int usedyads;      /*Is a free-dyad restriction in force?*/
  RLEBDM1D dyads;    /*...and if so, the free dyad set*/
  /*Edge list snapshot, rebuilt once per event for the degree families so
    that enumeration is over a stable array rather than a live tree.*/
  Vertex *eltail, *elhead;
  Edge nedges, elalloc;
} MoveSet;

/*An enumeration cursor over a move set.*/
typedef struct {
  const MoveSet *ms;
  Network *nwp;
  Dyad d;          /*MS_FREE/MS_DYAD: current dyad index*/
  RLERun hint;     /*MS_DYAD: run hint for NextRLEBDM1D*/
  Edge ei, ej;     /*Degree families: indices into the edge snapshot*/
  Vertex k;        /*Rewire families: the third node*/
  int variant;     /*MS_DEGREES: which of the two re-pairings*/
  int done;
} MoveIter;

MoveSet *MoveSetInit(int family, Network *nwp, double *rledyads);
void MoveSetDestroy(MoveSet *ms);
/*Refresh per-event caches (the edge snapshot).  Call after every accepted
  event, and before enumerating or sampling.*/
void MoveSetRefresh(MoveSet *ms, Network *nwp);

/*Size of the sampling superset: the number of candidates MoveSetSample
  draws from, uniformly.  Depends on the state only through the edge count,
  which the degree constraints hold fixed.*/
double MoveSetSupersetSize(const MoveSet *ms, Network *nwp);

/*Draw one candidate uniformly from the superset.  Returns TRUE if the
  candidate is a legal move (and fills mv), FALSE if it should be rejected.*/
Rboolean MoveSetSample(const MoveSet *ms, Network *nwp, Move *mv);

/*Enumeration.  Usage:
    MoveIter it; MoveIterInit(&it, ms, nwp);
    while(MoveIterNext(&it, &mv)){ ... }
*/
void MoveIterInit(MoveIter *it, const MoveSet *ms, Network *nwp);
Rboolean MoveIterNext(MoveIter *it, Move *mv);

#endif
