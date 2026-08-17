######################################################################
#
# constraints.R
#
# Licensed under the GNU General Public License version 3 or later
#
# Part of the R/ergmgp package
#
# Support for simulating EGPs on constrained sample spaces.
#
# An unconstrained EGP moves by toggling one dyad at a time.  Under a
# constraint such as ~degrees, *no* single toggle is legal (every toggle
# changes a degree), so the event model itself has to change: a move
# becomes a *set* of simultaneous toggles chosen so that the constrained
# statistic is preserved.  This file works out, for a given constraint
# formula and process, which move set to use and whether the combination
# is supported at all -- and says so explicitly rather than failing
# silently or quietly ignoring the constraint.
#
######################################################################


#Move-set family codes.  These must agree with the MS_* constants in
#src/moveset.h.
EGP_MS_FREE     <- 0L   #1 toggle, all dyads
EGP_MS_DYAD     <- 1L   #1 toggle, restricted dyad set (via RLE bit matrix)
EGP_MS_ODEGREES <- 2L   #2 toggles sharing a tail (out-degree preserving)
EGP_MS_IDEGREES <- 3L   #2 toggles sharing a head (in-degree preserving)
EGP_MS_DEGREES  <- 4L   #4 toggles (tetrad; degree preserving, undirected)
EGP_MS_EDGES    <- 5L   #2 toggles (one edge out, one non-edge in; edge count
                        #preserving, degree distribution free)

#Codes at or above this one are multi-toggle families, for which one event
#changes more than one edge variable.  Must agree with MS_MULTITOG in
#src/moveset.h.
EGP_MS_MULTITOG <- EGP_MS_ODEGREES


#Which constraints does ergmgp know how to turn into a move set?
#
#These are the count-preserving constraints: each fixes some function of the
#graph that no single toggle can leave alone, so each dictates its own move
#shape.  At most one of them can be in force at a time.
#
#`directed` says which networks the constraint applies to (NA = either).
#
#Dyad-level constraints are not listed here: they are detected
#structurally (see EGP_constraints below), so any present or future ergm
#constraint that exposes a free_dyads element works automatically.
EGP_moveset_constraints <- list(
  odegrees   = list(family=EGP_MS_ODEGREES, directed=TRUE,  ntoggles=2L),
  idegrees   = list(family=EGP_MS_IDEGREES, directed=TRUE,  ntoggles=2L),
  degrees    = list(family=EGP_MS_DEGREES,  directed=FALSE, ntoggles=4L),
  nodedegrees= list(family=EGP_MS_DEGREES,  directed=FALSE, ntoggles=4L),
  edges      = list(family=EGP_MS_EDGES,    directed=NA,    ntoggles=2L)
)


#Process capability table.
#
#  thin       - is the transition rate bounded above by a constant, so that
#               the process can be simulated exactly by uniformization
#               ("thinning")?  See EGPConstraintSupport() for the bounds.
#  multitog   - is the process well defined when a single event toggles more
#               than one dyad?  The continuum STERGMs are not: a rewire is
#               simultaneously a formation and a dissolution, so the
#               potential difference has no principled split into formation
#               and dissolution parts and separability breaks down.
EGP_process_caps <- list(
  LERGM     = list(thin=TRUE,  multitog=TRUE,  bound="A"),
  CRSAOM    = list(thin=FALSE, multitog=TRUE,  bound=NA_character_),
  CI        = list(thin=TRUE,  multitog=TRUE,  bound="A"),
  DS        = list(thin=TRUE,  multitog=TRUE,  bound="A/(|H|exp(pot))"),
  CDCSTERGM = list(thin=FALSE, multitog=FALSE, bound=NA_character_),
  CFCSTERGM = list(thin=FALSE, multitog=FALSE, bound=NA_character_),
  CSTERGM   = list(thin=FALSE, multitog=FALSE, bound=NA_character_),
  CTERGM    = list(thin=FALSE, multitog=TRUE,  bound=NA_character_)
)


#Report which constraints each EGP supports, and by which engine.  See
#man/EGPConstraintSupport.Rd; this package does not use roxygen.
EGPConstraintSupport<-function(process=NULL){
  cons<-c("~. (unconstrained)","dyad-level","~edges","~odegrees","~idegrees","~degrees")
  fams<-c(EGP_MS_FREE,EGP_MS_DYAD,EGP_MS_EDGES,EGP_MS_ODEGREES,EGP_MS_IDEGREES,EGP_MS_DEGREES)
  procs<-if(is.null(process)) names(EGP_process_caps) else match.arg(process,names(EGP_process_caps))
  out<-do.call(rbind,lapply(procs,function(p){
    cap<-EGP_process_caps[[p]]
    do.call(rbind,lapply(seq_along(cons),function(i){
      multi<-fams[i]>=EGP_MS_MULTITOG
      supported<-(!multi)||cap$multitog
      #DS is a special case: its rate involves the size of the move set, which
      #for a multi-toggle family can only be had by enumerating it, so thinning
      #is unavailable there however bounded the rate is.
      thin<-cap$thin && !(p=="DS" && multi)
      data.frame(process=p, constraint=cons[i],
                 supported=supported,
                 engine=if(!supported) "-" else if(thin) "thinning (exact)" else "enumeration",
                 note=if(!supported) "multi-toggle moves break formation/dissolution separability"
                      else if(!multi) ""
                      else if(thin) paste0("rate bound = ",cap$bound)
                      else if(p=="DS") "rate needs the move count, so it must be enumerated (cheaply)"
                      else "no rate bound; enumerated (slow for multi-toggle moves)",
                 stringsAsFactors=FALSE)
    }))
  }))
  rownames(out)<-NULL
  print(out,right=FALSE)
  cat("\nDyad-level constraints are those exposing a free-dyad set, e.g. blocks,\n",
      "observed, fixedas, fixallbut, egocentric, and combinations thereof.\n",
      "They compose with the count-preserving constraints (intersection).\n",sep="")
  invisible(out)
}


#Refuse to silently ignore unrecognised arguments.
#
#simEGP() historically ended in `...`, documented as "not currently used".
#That made every typo and every unsupported option a silent no-op: passing
#constraint=~odegrees (singular) ran the simulation unconstrained and gave
#no indication that anything had been ignored.  This is the single most
#dangerous failure mode in the package, so we now stop instead.
EGP_check_dots<-function(...){
  extra<-list(...)
  if(length(extra)==0L) return(invisible(NULL))
  nm<-names(extra)
  if(is.null(nm)) nm<-rep("",length(extra))
  hint<-character(0)
  if("constraint"%in%nm)
    hint<-c(hint,"  Did you mean `constraints` (plural)?  ergmgp follows ergm's spelling.")
  if(any(nm==""))
    hint<-c(hint,"  Unnamed extra arguments are never used; check for a misplaced comma.")
  stop("Unused argument(s) passed to simEGP(): ",
       paste(ifelse(nm=="","<unnamed>",nm),collapse=", "),".\n",
       if(length(hint)) paste(hint,collapse="\n") else
         "  These were silently ignored by earlier versions; they are now an error.",
       call.=FALSE)
}


#Parse and validate a constraint formula for a given process.
#
#Arguments:
#  constraints - a one-sided constraint formula (as for ergm), or NULL/~.
#  nw          - the network being simulated
#  proc        - the process name (character)
#
#Return value:
#  A list with elements
#    family   - move-set family code (EGP_MS_*)
#    rlebdm   - free-dyad RLE bit matrix as a numeric vector for C, or NULL
#    ntoggles - number of toggles per move
#    label    - human-readable description, for messages
#
#This function is deliberately loud.  Every unsupported combination raises
#an error naming the process, the constraint, and the reason, because the
#failure mode we are guarding against is a simulation that silently runs
#unconstrained and produces plausible-looking but wrong results.
EGP_constraints<-function(constraints, nw, proc){
  cap<-EGP_process_caps[[proc]]
  #No constraint supplied: the classical single-toggle EGP
  if(is.null(constraints))
    return(list(family=EGP_MS_FREE, rlebdm=NULL, ntoggles=1L, label="~. (unconstrained)"))
  if(!inherits(constraints,"formula"))
    stop("`constraints` must be a one-sided formula, e.g. constraints=~degrees.  ",
         "You supplied an object of class ",sQuote(class(constraints)[1]),".")
  if(length(constraints)==3L)
    stop("`constraints` must be one-sided (constraints=~degrees), not two-sided.")

  #Parse via ergm, so we accept exactly what ergm accepts
  conlist<-ergm::ergm_conlist(constraints, nw)
  cnames<-setdiff(names(conlist),c(".","")) #"." is the implicit no-op constraint
  if(length(cnames)==0L)
    return(list(family=EGP_MS_FREE, rlebdm=NULL, ntoggles=1L, label="~. (unconstrained)"))

  #Split into count-preserving vs dyad-level.  Dyad-level constraints are
  #detected structurally: they are exactly those carrying a free_dyads
  #element, which is how ergm itself distinguishes them.
  #
  #Note that ergm has already resolved implication among the count-preserving
  #constraints for us: ~degrees implies ~edges, so ergm_conlist(~edges+degrees)
  #returns only `degrees` and we never see the redundant pair.
  cntcon<-intersect(cnames,names(EGP_moveset_constraints))
  othercon<-setdiff(cnames,cntcon)
  isdyadlevel<-vapply(othercon,function(k) !is.null(conlist[[k]]$free_dyads), logical(1))
  unsupported<-othercon[!isdyadlevel]
  if(length(unsupported)>0L)
    stop("ergmgp cannot simulate under the constraint(s) ",
         paste0("~",unsupported,collapse=", "),".\n",
         "  Supported: dyad-level constraints (blocks, observed, fixedas, fixallbut,\n",
         "  egocentric, ...) and the count-preserving constraints ",
         paste0("~",names(EGP_moveset_constraints),collapse=", "),".\n",
         "  Call EGPConstraintSupport() for the full table.")
  if(length(cntcon)>1L)
    stop("Only one count-preserving constraint can be imposed at a time; you gave ",
         paste0("~",cntcon,collapse=", "),".\n",
         "  Each of them fixes a different graph statistic and so dictates a different\n",
         "  move shape; ergmgp has no move set preserving two of them at once.")

  #Count-preserving constraints: check directedness and process compatibility
  family<-EGP_MS_FREE; ntog<-1L; label<-character(0)
  if(length(cntcon)==1L){
    spec<-EGP_moveset_constraints[[cntcon]]
    #directed=NA means the constraint applies to both kinds of network
    if(!is.na(spec$directed)){
      if(spec$directed && !is.directed(nw))
        stop("Constraint ~",cntcon," applies to directed networks only, but this network ",
             "is undirected.  For undirected networks use constraints=~degrees, which ",
             "preserves the full degree sequence.")
      if(!spec$directed && is.directed(nw))
        stop("Constraint ~",cntcon," applies to undirected networks only, but this network ",
             "is directed.  For directed networks use ~odegrees or ~idegrees.")
    }
    if(!cap$multitog)
      stop("Process ",proc," cannot be simulated under ~",cntcon,".\n",
           "  ~",cntcon," requires each event to toggle ",spec$ntoggles," dyads at once, but ",
           proc," is a\n  continuum STERGM: such a move is simultaneously a formation and a ",
           "dissolution, so the\n  potential difference has no principled split into formation and ",
           "dissolution parts\n  and separability breaks down.  Use LERGM, CI, DS, CRSAOM or CTERGM, ",
           "or restrict\n  yourself to dyad-level constraints.  See EGPConstraintSupport().")
    family<-spec$family; ntog<-spec$ntoggles; label<-paste0("~",cntcon)
  }

  #Dyad-level constraints: extract the free dyad set for C.
  #
  #ergm always adds an implicit ".attributes" constraint whose free-dyad set is
  #just the dyads the network has at all (e.g. the upper triangle when
  #undirected).  That is not a restriction -- the move sets already respect it
  #-- so we only carry a dyad set through to C when it is a proper subset of
  #the trivially free dyads.  Anything genuinely restrictive is kept.
  rle<-NULL
  if(length(othercon)>0L){
    fd<-ergm::as.rlebdm(conlist)
    if(is.null(fd))
      stop("Could not derive a free-dyad set from constraints ",
           paste0("~",othercon,collapse=", "),".")
    nfree<-sum(as.numeric(fd$lengths[fd$values]))
    nn<-network.size(nw)
    ntriv<-if(is.directed(nw)) nn*(nn-1) else nn*(nn-1)/2
    if(nfree==0)
      stop("The constraint ",paste0("~",setdiff(othercon,".attributes"),collapse=", "),
           " leaves no free dyads at all, so no event can ever occur.\n",
           "  Check the constraint's arguments: this is almost always a specification\n",
           "  error rather than an intended frozen network.  (For `blocks`, note that\n",
           "  `levels2` selects the block pairs to *forbid*.)")
    if(nfree<ntriv){
      rle<-as.double(ergm::to_ergm_Cdouble(fd))
      if(family==EGP_MS_FREE){ family<-EGP_MS_DYAD; ntog<-1L }
      label<-c(label,paste0("~",setdiff(othercon,".attributes"),collapse=", "))
    }
  }
  if(length(label)==0L) label<-"~. (unconstrained)"

  list(family=family, rlebdm=rle, ntoggles=ntog,
       label=paste(label,collapse=" + "))
}


#Decide which simulation engine to use, and say so.
#
#Arguments:
#  con     - the object returned by EGP_constraints()
#  proc    - process name
#  engine  - user request: "auto", "thinning", or "enumeration"
#  verbose - logical; report the choice?
#
#Return value: 0L for enumeration, 1L for thinning.
EGP_engine<-function(con, proc, engine=c("auto","thinning","enumeration"), verbose=FALSE){
  engine<-match.arg(engine)
  cap<-EGP_process_caps[[proc]]
  if(engine=="thinning" && !cap$thin)
    stop("engine=\"thinning\" is not available for process ",proc,".\n",
         "  Thinning (uniformization) requires the transition rate to be bounded above by a\n",
         "  constant, which holds for LERGM, CI and DS but not for ",proc,".  Use\n",
         "  engine=\"enumeration\" (exact, but slow for a multi-toggle move set) or a\n",
         "  different process.")
  #The differential stability rate is A/(|H| exp(pot)), where |H| is the number
  #of legal moves.  For the dyad families |H| is known in closed form, but for
  #the multi-toggle families counting the legal moves is exactly the enumeration
  #we were trying to avoid, so thinning buys nothing and cannot be made exact.
  if(proc=="DS" && con$family>=EGP_MS_MULTITOG){
    if(engine=="thinning")
      stop("engine=\"thinning\" is not available for DS under ",con$label,".\n",
           "  The DS rate is A/(|H| exp(pot)), and |H| is the number of legal moves, which\n",
           "  for a multi-toggle constraint can only be obtained by enumerating them.\n",
           "  Use engine=\"enumeration\" (DS enumerates cheaply: because its rate does not\n",
           "  depend on the move, no change statistics are computed during the sweep).")
    engine<-"enumeration"
  }
  #"auto" keeps the unconstrained process on enumeration, which is what every
  #released version of this package has done: switching its engine silently
  #would change results people already have, for a speedup they did not ask
  #for.  Thinning is available there via engine="thinning", and is typically
  #one to two orders of magnitude faster.  For the constrained families there
  #is no prior behaviour to preserve, and enumeration is prohibitively slow, so
  #thinning is the default wherever the process admits a rate bound.
  use.thin<-switch(engine,
                   auto=cap$thin && con$family!=EGP_MS_FREE,
                   thinning=TRUE,
                   enumeration=FALSE)
  if(verbose){
    message("EGP: process ",proc,", constraint ",con$label,", ",
            con$ntoggles," toggle(s) per event, engine ",
            if(use.thin) "thinning (exact uniformization)" else "enumeration",".")
    if(!use.thin && con$ntoggles>1L)
      message("EGP: note - enumerating a ",con$ntoggles,
              "-toggle move set is expensive; consider a process that supports thinning ",
              "(LERGM, CI, DS).  See EGPConstraintSupport().")
  }
  as.integer(use.thin)
}
