/*
######################################################################
#
# simEGP.c
#
# copyright (c) 2023, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/4/23
# Licensed under the GNU General Public License version 3 or later.
#
# Part of the R/ergmgp package
#
# This file contains functions for simulating realizations from ERGM
# generating processes.
#
######################################################################
*/

#include "egp.h"
#include "ergm_model.h"
#include "ergm_state.h"
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_util.h"
#include <R_ext/Utils.h>
#include "utils.h"
#include "moveset.h"

/*
  Compute the change statistics, relative potential, and log transition rate
  for a single move.

  A move is a set of dyads to be toggled simultaneously.  For the ordinary
  one-dyad move we call the term c_functions directly, which is cheap because
  it never touches the network.  For multi-dyad moves (needed by the degree
  preserving constraints) we have to go through ergm's ChangeStats(), which
  applies the toggles one at a time and then rolls them back: terms such as
  gwesp supply only a c_function, so there is no other way to get a joint
  change statistic out of them.  That rollback is why multi-toggle moves cost
  substantially more than the toggle count alone would suggest.

  Returns the log rate, and writes the relative potential(s) to relpot.
*/
static double score_move(const Move *mv, int egp, Model *m, Network *nwp,
                         double *coef, int cooffset, double lcrate,
                         double lHneigh, const double *pot, double *relpot,
                         double *newRow){
  int j;
  relpot[0] = relpot[1] = 0.0;
  memset(newRow, 0.0, m->n_stats*sizeof(double));
  if(mv->ntoggles == 1){
    Vertex t = mv->tails[0], h = mv->heads[0];
    EXEC_THROUGH_TERMS_INTO(m, newRow, {
      if(mtp->c_func){
        ZERO_ALL_CHANGESTATS();
        (*(mtp->c_func))(t, h, mtp, nwp, IS_OUTEDGE(t,h));
      }else if(mtp->d_func){
        (*(mtp->d_func))(1, &t, &h, mtp, nwp);
      }
      addonto(dstats, mtp->dstats, N_CHANGE_STATS);
    });
  }else{
    ChangeStats(mv->ntoggles, (Vertex *)mv->tails, (Vertex *)mv->heads, nwp, m);
    memcpy(newRow, m->workspace, m->n_stats*sizeof(double));
  }
  /*Compute relative potential (or potentials, for CSTERGMs)*/
  if(egp==EGP_CDCSTERGM){        /*Constant dissolution CSTERGMs*/
    for(j=0;j<m->n_stats;j++)
      relpot[0]+=coef[j+1]*newRow[j];
    relpot[1]=coef[0];
  }else if(egp==EGP_CFCSTERGM){  /*Constant formation CSTERGMs*/
    for(j=0;j<m->n_stats;j++)
      relpot[1]+=coef[j+1]*newRow[j];
    relpot[0]=coef[0];
  }else if(egp==EGP_CSTERGM){     /*For CSTERGMs, consider form and diss*/
    for(j=0;j<cooffset;j++)
      relpot[0]+=coef[j]*newRow[j];
    for(j=cooffset;j<m->n_stats;j++)
      relpot[1]+=coef[j]*newRow[j];
  }else{                          /*Standard calculation for most models*/
    for(j=0;j<m->n_stats;j++)
      relpot[0]+=coef[j]*newRow[j];
  }
  /*Turn the potential difference into a rate, per the process definition*/
  switch(egp){
    /*Longitudinal ERGM: log rate = log(A) - log1p(-potdiff)*/
    case EGP_LERGM:      return lcrate - logspace_add(0.0,-relpot[0]);
    /*Competing rate SAOM: log rate = log(A) + targpot*/
    case EGP_CRSAOM:     return lcrate + relpot[0] + pot[0];
    /*Change inhibition process: log rate = log(A) + min(0, potdiff)*/
    case EGP_CI:         return lcrate + MIN(0.0,relpot[0]);
    /*Differential stability: log rate = log(A) - log |H| - curpot*/
    case EGP_DS:         return lcrate - lHneigh - pot[0];
    /*Continuum STERGMs: log rate = log(A) + (edge ? potdiff_d : potdiff_f)*/
    case EGP_CDCSTERGM:
    case EGP_CFCSTERGM:
    case EGP_CSTERGM:    return lcrate + (mv->isedge ? relpot[1] : relpot[0]);
    /*Continuum TERGM: log rate = log(A) + potdiff*/
    case EGP_CTERGM:     return lcrate + relpot[0];
  }
  return R_NegInf;
}

/*
  The per-move rate bound used by the thinning (uniformization) engine.

  Thinning is exact whenever the transition rate of every move is bounded
  above by a constant Lambda: we subordinate the process to a Poisson process
  of rate N*Lambda over a fixed candidate superset of size N, and accept each
  candidate with probability rate/Lambda.  Rejections are self-loops, which
  advance time without changing state.  Returns R_PosInf when no bound
  exists, in which case the caller must enumerate instead.
*/
static double rate_bound(int egp, double lcrate, double lHneigh, const double *pot){
  switch(egp){
    case EGP_LERGM: return lcrate;                       /*A/(1+exp(-d)) <= A*/
    case EGP_CI:    return lcrate;                       /*A*min(1,exp(d)) <= A*/
    case EGP_DS:    return lcrate - lHneigh - pot[0];    /*Rate is move-independent*/
    default:        return R_PosInf;
  }
}


/*
  Simulate a single trajectory from an ERGM generating process (EGP); several
  different types of EGPs are supported.  This function simulates a trajectory
  of specified length and returns the final state (along with, optionally, history
  information); length may be specified by the number of events drawn (time will
  be random) or by the time taken (event count will be random).  It is important
  to note that event updates are *not* random times, and hence trajectories
  terminated by event count will not be in equilibrium.

  Arguments:
    segp - an ERGM generating process type code, saying which process should be used
           for the simulation
    mstate - an input ergm state, from ergm_state(), containing the model, graph state,
             and anything else needed
    scoef - coefficient vector for the model; for constant rate CSTERGMs, the first
            element should be the constant rate (formation or dissolution).  For
            general CSTERGMs, this should be a vector of formation coefficients
            followed by the dissolution coefficients
    scooffset - coefficient offsets for CSTERGMs.  For constant rate CSTERGMs, this
                should be 1, and the first element should be whichever rate is
                constant; for general CSTERGMs, this should be the number of formation
                parameters
    stmax - maximum simulation time; if constraining by events, this should be negative
            (meaning unlimited or infinite)
    sevmax - maximum number of events to simulate; if constraining by time, this should be
               negative (meaning unlimited or infinite)
    slcrate - log "collision rate" or intrinsic rate coefficient; this is the maximum
              transition rate as the log potential difference goes to infinity for some
              processes, and for others simply sets the timescale of dynamics
    spot - "negative effective energy" (aka ERGM potential) at onset; this allows
              the potential to be tracked
    slogtime - logical; track time on log scale?  This is purely internal; using it
               can potentially avoid overflow/underflow errors for models with extreme
               rates, but also increases computational burden somewhat
    slasttog - logical; or every edge variable, track the time of its last toggle
               (i.e., from off to on, or on to off)?  If TRUE, an n x n matrix is
               included in the return value conaining the last toggle time
    slttoff - n x n matrix of initial last toggle times, in case we are resuming a
               simulation already in progress
    skeephist - logical; should we keep and return the entire event history?  (Not for
                the faint of heart or those short of RAM.)
    sverbose - if positive, we will return progress messages; events are shown every
               verbose events (so e.g., verbose=100 is a good idea, and verbose=1 is not)
    sfamily - move-set family code (see moveset.h); 0 for the classical unconstrained
              single-toggle process
    srledyads - free dyad set for dyad-level constraints, in to_ergm_Cdouble form, or
                NULL if there is no dyad-level restriction
    sthin - logical; use the thinning (uniformization) engine rather than enumerating
            the move set?  Only valid for processes with a bounded rate

  Return value:
    A named list with elements
      state: the final state, in ergm_state() form
      evcount: the final event or step count (number of observed transitions)
      simtime: the final simulation time (which may not be the time of the last event,
               depending on whether events or total time were left random)
      potential: final ERGM potential
      lasttog: if used, an n x n matrix containing the last toggle time for each edge
               variable
      evhist: if used, the complete event history (as an event by 4 matrix, with cols
               sender, receiver, time, and formation indicator (1=onset, 0=terminus)
      proposals: number of candidate moves drawn (thinning engine only)
*/
SEXP simEGP_R(SEXP segp, SEXP mstate, SEXP scoef, SEXP scooffset, SEXP stmax, SEXP sevmax, SEXP slcrate, SEXP spot, SEXP slogtime, SEXP slasttog, SEXP slttoff, SEXP skeephist, SEXP sverbose, SEXP sfamily, SEXP srledyads, SEXP sthin){
  int pc=0,verbose,i,j,egp,cooffset,logtime,lasttog,keephist,family,thin;
  double *coef,tmax,evmax,lcrate,ecount,simtime,pot[2],lHneigh,*ltt,*dp,*rledyads;
  double tm,relpot[2],ltotrat=R_NegInf,selpot[2];
  double ncand=0.0,propcount=0.0;
  double *pehist,*pdp;
  element *ehist;
  SEXP outl,stats,sltt,sehist;
  R_xlen_t nev,xi;
  Move mv,selmv;
  selpot[0]=selpot[1]=R_NegInf;
  ehist=NULL;
  sltt=NULL;
  ltt=NULL;
  memset(&selmv,0,sizeof(Move));

  /*Ensure that inputs are processed correctly*/
  PROTECT(segp=coerceVector(segp,INTSXP));           /*EGP type*/
  egp=INTEGER(segp)[0];
  UNPROTECT(1);
  PROTECT(scoef=coerceVector(scoef,REALSXP)); pc++;  /*Model coefficients*/
  coef=REAL(scoef);
  PROTECT(stmax=coerceVector(stmax,REALSXP));        /*Max simulation time*/
  tmax=REAL(stmax)[0];
  UNPROTECT(1);
  PROTECT(sevmax=coerceVector(sevmax,REALSXP));      /*Max events*/
  evmax=REAL(sevmax)[0];
  UNPROTECT(1);
  PROTECT(slcrate=coerceVector(slcrate,REALSXP));    /*Log collision rate*/
  lcrate=REAL(slcrate)[0];
  UNPROTECT(1);
  PROTECT(spot=coerceVector(spot,REALSXP));          /*ERGM potential*/
  pot[0]=REAL(spot)[0];                                /*Either total or formation pot*/
  pot[1]=REAL(spot)[1];                                /*Either unused or dissolution pot*/
  UNPROTECT(1);
  PROTECT(sverbose=coerceVector(sverbose,INTSXP));   /*Verbosity*/
  verbose=INTEGER(sverbose)[0];
  UNPROTECT(1);
  PROTECT(slogtime=coerceVector(slogtime,INTSXP));   /*Use log time*/
  logtime=INTEGER(slogtime)[0];
  UNPROTECT(1);
  PROTECT(slasttog=coerceVector(slasttog,INTSXP));   /*Track last toggles*/
  lasttog=INTEGER(slasttog)[0];
  UNPROTECT(1);
  PROTECT(skeephist=coerceVector(skeephist,INTSXP)); /*Keep event history*/
  keephist=INTEGER(skeephist)[0];
  UNPROTECT(1);
  PROTECT(scooffset=coerceVector(scooffset,INTSXP));  /*CSTERGM form coef offset*/
  cooffset=INTEGER(scooffset)[0];
  UNPROTECT(1);
  PROTECT(sfamily=coerceVector(sfamily,INTSXP));     /*Move-set family*/
  family=INTEGER(sfamily)[0];
  UNPROTECT(1);
  PROTECT(sthin=coerceVector(sthin,INTSXP));         /*Thinning engine?*/
  thin=INTEGER(sthin)[0];
  UNPROTECT(1);
  if(srledyads==R_NilValue){
    rledyads=NULL;
  }else{
    PROTECT(srledyads=coerceVector(srledyads,REALSXP)); pc++;
    rledyads=REAL(srledyads);
  }

  /*Initialize other key things*/
  GetRNGstate();
  ErgmState *s=ErgmStateInit(mstate, ERGM_STATE_NO_INIT_PROP);  /*Set up the ERGM state*/
  Model *m = s->m;              /*Extract the model from the state*/
  Network *nwp=s->nwp;           /*Extract the network from the state*/
  int directed_flag = nwp->directed_flag;  /*Is this directed?*/
  Vertex n_nodes = nwp->nnodes;
  if(tmax<0.0)                 /*Deal w/R's prohibition on passing Inf values*/
    tmax=R_PosInf;
  if(evmax<0.0)
    evmax=R_PosInf;
  if((egp==EGP_CDCSTERGM)||(egp==EGP_CFCSTERGM)) /*Set coef offset for constant rate CSTERGMs*/
    cooffset=1;
  /*Set up the move set*/
  MoveSet *ms=MoveSetInit(family, nwp, rledyads);
  /*Log Hamming neighbourhood size.  For the differential stability process this
    is part of the rate, so it must reflect the constrained move set rather than
    the full dyad space.*/
  if(family==MS_DYAD)
    lHneigh=log((double)ms->dyads.ndyads);
  else
    lHneigh=log(n_nodes*(n_nodes-1.0)/(2.0-directed_flag));
  if(lasttog){     /*If storing last toggle times, initialize the n x n timing matrix*/
    if(slttoff==R_NilValue){  /*Start from scratch*/
      PROTECT(sltt=allocMatrix(REALSXP,n_nodes,n_nodes)); pc++;
      ltt=REAL(sltt);
      for(i=0;i<n_nodes;i++)
        for(j=0;j<n_nodes;j++)
          ltt[i+j*n_nodes]=0.0;
    }else{                      /*Use the one we were given*/
      PROTECT(sltt=coerceVector(slttoff,REALSXP)); pc++;
      if(LENGTH(sltt)!=n_nodes*n_nodes){
        error("We were given a last toggle time matrix of total length %ld, but it should have been %ld.  Stopping.\n", (long)LENGTH(sltt), (long)(n_nodes*n_nodes));
      }
      ltt=REAL(sltt);
    }
  }

  /*Memory to store one row of changescores*/
  SEXP snewrow;
  PROTECT(snewrow=allocVector(REALSXP,m->n_stats)); pc++;
  double *newRow = REAL(snewrow);
  memset(newRow,0.0,m->n_stats*sizeof(double));

  /*Memory to store a matrix of log rates.  Only the single-toggle families use
    this: the degree-preserving move sets can have O(E^2) members, far too many
    to hold, so they select by the Gumbel-max trick instead.*/
  SEXP slrates;
  int storerates = (family==MS_FREE)||(family==MS_DYAD);
  PROTECT(slrates=allocVector(REALSXP, storerates ? n_nodes*n_nodes : 1)); pc++;
  double *lrates = REAL(slrates);
  if(storerates)
    memset(lrates,0.0,n_nodes*n_nodes*sizeof(double));

  /*Memory for ergm stat state (needs to be protected)*/
  PROTECT(stats=allocVector(REALSXP,m->n_stats)); pc++;
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));

  /*Run the event simulation until time expires or we run out of events*/
  if(logtime){
    tm=R_NegInf;   /*Store time on log scale, to deal with extreme rate cases*/
    if(tmax<R_PosInf)
      tmax=log(tmax);
    simtime=R_NegInf;
  }else{
    tm=0.0;        /*...or just use standard scale, which is faster*/
    simtime=0.0;
  }
  ecount=0.0;
  if(verbose){
    Rprintf("Initializing simulation: max events=%.0f, max time=%0f, initial pot=(%f,%f)\n", evmax, tmax, pot[0],pot[1]);
  }
  while((ecount<evmax)&&(tm<tmax)){
    if((verbose)&&(((int)(ecount))%verbose==0)){
      Rprintf("event=%.0f, t=%f, pot=(%f,%f)\n",ecount,tm,pot[0],pot[1]);
    }
    R_CheckUserInterrupt();
    int havemove=0;

    if(thin){
      /*----------------------------------------------------------------
        Thinning (uniformization) engine.

        We draw candidates uniformly from a fixed superset of the move
        set whose size N does not depend on the current state, wait
        Exp(N*Lambda) between candidates, and accept each with
        probability rate/Lambda.  Illegal or rejected candidates are
        self-loops: time advances, the state does not.  This is exact,
        and costs one change-statistic evaluation per candidate rather
        than one per member of the move set.
        ----------------------------------------------------------------*/
      double lLambda=rate_bound(egp, lcrate, lHneigh, pot);
      ncand=MoveSetSupersetSize(ms, nwp);
      if(!R_FINITE(lLambda)||ncand<=0.0){
        error("Internal error: the thinning engine was invoked for a process with no rate bound, or an empty move set.  Please report this.\n");
      }
      ltotrat=log(ncand)+lLambda;
      /*Advance time to the next candidate*/
      if(logtime)
        tm=logspace_add(tm, log(-log(runif(0.0,1.0)))-ltotrat);
      else
        tm+=rexp(exp(-ltotrat));
      if(tm>=tmax) break;
      propcount++;
      if(MoveSetSample(ms, nwp, &mv)){
        double lrate=score_move(&mv, egp, m, nwp, coef, cooffset, lcrate, lHneigh, pot, relpot, newRow);
        /*Accept with probability rate/Lambda*/
        if(log(unif_rand())<lrate-lLambda){
          selmv=mv;
          selpot[0]=relpot[0];
          selpot[1]=relpot[1];
          havemove=1;
        }
      }
      if(!havemove)
        continue;   /*Self-loop: time has advanced, state has not*/
    }else{
      /*----------------------------------------------------------------
        Enumeration engine.

        Walk the entire legal move set, accumulate the total rate, and
        select the winner.  Exact for every process, but the move set is
        O(E^2) for ~degrees, which is why the thinning engine exists.
        ----------------------------------------------------------------*/
      MoveIter it;
      MoveSetRefresh(ms, nwp);
      MoveIterInit(&it, ms, nwp);
      ltotrat=R_NegInf;
      if(egp==EGP_DS){
        /*The differential stability rate does not depend on the move, so we
          need only count the legal moves and pick one uniformly (reservoir
          sampling).  Skipping the change statistics here is a large saving.*/
        double nmoves=0.0;
        while(MoveIterNext(&it, &mv)){
          nmoves++;
          if(unif_rand()*nmoves<1.0){ selmv=mv; havemove=1; }
        }
        if(nmoves<=0.0) break;   /*No legal move: the process is absorbed*/
        lHneigh=log(nmoves);
        ltotrat=lcrate-pot[0];   /*= nmoves * A/(|H| exp(pot)) with |H|=nmoves*/
        if(havemove)
          score_move(&selmv, egp, m, nwp, coef, cooffset, lcrate, lHneigh, pot, selpot, newRow);
      }else{
        double selmax=R_NegInf,lmax=R_NegInf;
        double nmoves=0.0;
        if(storerates){
          /*Small, bounded move set: cache the rates so that we can select by
            walking a cumulative sum, rather than paying two logarithms per
            move for a Gumbel draw.*/
          while(MoveIterNext(&it, &mv)){
            double lrate=score_move(&mv, egp, m, nwp, coef, cooffset, lcrate, lHneigh, pot, relpot, newRow);
            lrates[mv.tails[0]-1+(mv.heads[0]-1)*n_nodes]=lrate;
            if(lrate>lmax) lmax=lrate;
            nmoves++;
          }
          if(nmoves<=0.0) break;
          /*Total rate, computed with an explicit max shift for stability*/
          double sum=0.0;
          MoveIterInit(&it, ms, nwp);
          while(MoveIterNext(&it, &mv))
            sum+=exp(lrates[mv.tails[0]-1+(mv.heads[0]-1)*n_nodes]-lmax);
          ltotrat=lmax+log(sum);
          /*Select by walking the cumulative sum*/
          double u=unif_rand()*sum,acc=0.0;
          MoveIterInit(&it, ms, nwp);
          while(MoveIterNext(&it, &mv)){
            acc+=exp(lrates[mv.tails[0]-1+(mv.heads[0]-1)*n_nodes]-lmax);
            if(acc>=u){ selmv=mv; havemove=1; break; }
          }
          if(!havemove){ selmv=mv; havemove=1; }  /*Rounding guard*/
          /*Recover the winner's potential difference (one extra evaluation)*/
          score_move(&selmv, egp, m, nwp, coef, cooffset, lcrate, lHneigh, pot, selpot, newRow);
        }else{
          /*Large move set: select by the Gumbel-max trick, which needs no
            storage.  We draw a Gumbel(0) deviate for each move and take the
            arg max of log(rate) + G, which can be done as we go.*/
          while(MoveIterNext(&it, &mv)){
            double lrate=score_move(&mv, egp, m, nwp, coef, cooffset, lcrate, lHneigh, pot, relpot, newRow);
            ltotrat=logspace_add(ltotrat,lrate);
            double selval=lrate-log(-log(runif(0.0,1.0)));
            if(selval>selmax){
              selmax=selval;
              selmv=mv;
              selpot[0]=relpot[0];
              selpot[1]=relpot[1];
              havemove=1;
            }
            nmoves++;
          }
          if(nmoves<=0.0) break;   /*No legal move: the process is absorbed*/
        }
      }
      /*Draw the next event time*/
      if(logtime)
        tm=logspace_add(tm, log(-log(runif(0.0,1.0)))-ltotrat); /*Log version*/
      else
        tm+=rexp(exp(-ltotrat)); /*(Non log version) Note that they use scale rather than rate here...*/
    }

    /*If the next event would occur before the end of the period, make the transition*/
    if((tm<tmax)&&havemove){
      ecount++;                   /*Increment the event count*/
      simtime=tm;                 /*Update the time of the last event*/
      pot[0] += selpot[0];        /*Formation (or general) potential*/
      pot[1] += selpot[1];        /*Dissolution potential*/
      for(i=0;i<selmv.ntoggles;i++){
        int wasedge=IS_OUTEDGE(selmv.tails[i],selmv.heads[i]);
        TOGGLE(selmv.tails[i],selmv.heads[i]);  /*Toggle the winning move!*/
        if(lasttog){                /*Store the update time, if desired*/
          ltt[selmv.tails[i]-1+(selmv.heads[i]-1)*n_nodes]=(logtime ? exp(tm) : tm);
        }
        if(keephist){
          dp=(double *)R_alloc(3,sizeof(double));
          dp[0]=(double)selmv.tails[i];
          dp[1]=(double)selmv.heads[i];
          dp[2]=1.0-(double)wasedge;
          ehist=enqueue(ehist,(logtime ? exp(tm) : tm), (void *)dp);
        }
      }
      if(!thin)
        MoveSetRefresh(ms, nwp);   /*Only enumeration reads the edge snapshot*/
    }
  }
  if(tm>tmax){  /*If we truncated at time tmax, update simulation time accordingly*/
    simtime=tmax;
  }
  if(logtime)  /*If we used logtime, transform back*/
    simtime=exp(simtime);

  /*Save our stuffs*/
  const char *outputnams[] = {"state","evcount","simtime","potential","lasttog","evhist","proposals",""};
  PROTECT(outl=mkNamed(VECSXP,outputnams)); pc++;
  s->stats=REAL(stats);
  SET_VECTOR_ELT(outl,0,ErgmStateRSave(s));
  SEXP secount,ssimtime,sprop;
  PROTECT(secount=allocVector(REALSXP,1)); pc++;
  PROTECT(ssimtime=allocVector(REALSXP,1)); pc++;
  PROTECT(sprop=allocVector(REALSXP,1)); pc++;
  REAL(secount)[0]=ecount;
  REAL(ssimtime)[0]=simtime;
  REAL(sprop)[0]=propcount;
  REAL(spot)[0]=pot[0];
  REAL(spot)[1]=pot[1];
  SET_VECTOR_ELT(outl,1,secount);
  SET_VECTOR_ELT(outl,2,ssimtime);
  SET_VECTOR_ELT(outl,3,spot);
  SET_VECTOR_ELT(outl,6,sprop);
  if(lasttog){                         /*Save the toggle matrix*/
    SET_VECTOR_ELT(outl,4,sltt);
  }
  if(keephist){                        /*Save the event history matrix*/
    /*One event can toggle several dyads, so count the queue rather than
      assuming one row per event.*/
    nev=0;
    {
      element *ep=ehist;
      while(ep!=NULL){ nev++; ep=ep->next; }
    }
    PROTECT(sehist=allocMatrix(REALSXP,nev,4)); pc++;
    pehist=REAL(sehist);
    for(xi=nev-1;xi>=0;xi--){
      pehist[xi]=ehist->val;                          /*Event time*/
      pdp=(double *)(ehist->dp);
      pehist[xi+nev]=pdp[0];                          /*Tail*/
      pehist[xi+2*nev]=pdp[1];                        /*Head*/
      pehist[xi+3*nev]=pdp[2];                        /*Onset/terminus indicator*/
      ehist=ehist->next;
    }
    SET_VECTOR_ELT(outl,5,sehist);
  }
  /*Clean up and return*/
  MoveSetDestroy(ms);
  ErgmStateDestroy(s);
  PutRNGstate();
  UNPROTECT(pc);
  return outl;
}
