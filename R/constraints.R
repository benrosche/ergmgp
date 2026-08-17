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


#Does the number of *legal* moves take the same value in every state the
#process can reach?
#
#This matters for one process only.  DS's rate is A/(|H| exp(pot)), so its total
#exit rate is A exp(-pot) and its jump chain is uniform over the |H| legal
#moves.  The move graph is symmetric, so that jump chain has stationary
#distribution nu(x) proportional to |H|(x), and for a CTMC pi = nu/q, giving
#
#    pi_DS(x)  proportional to  |H|(x) exp(pot(x))
#
#which is the intended ERGM if and only if |H| is constant.  Counting the moves
#each family enumerates (see MoveIterNext in src/moveset.c):
#
#  MS_FREE       all dyads                                   D
#  MS_DYAD       the free dyads                              |free|
#  MS_EDGES      (edge, non-edge) pairs                      E(D - E)
#  MS_ODEGREES   (edge, new head) pairs                      sum_t d_t(n - 1 - d_t)
#  MS_IDEGREES   (edge, new tail) pairs                      sum_h d_h(n - 1 - d_h)
#  MS_DEGREES    node-disjoint edge pairs whose two
#                re-pairings are both currently absent       state-dependent
#
#Every count but the last is fixed by exactly what its constraint holds fixed,
#so |H| cannot move.  This survives composition with a dyad-level constraint:
#frozen dyads never toggle, so the free edges and free non-edges at each node
#are themselves constant, and the products above just run over the free part.
#
#The tetrad is different in kind.  It is legal only if the two dyads it would
#*form* are currently absent, which is a fact about the configuration and not
#about the degree sequence -- so |H| genuinely varies (250-310 over equilibrium
#draws from a 14-node, 28-edge model), and DS under ~degrees equilibrates to a
#distribution tilted by |H| rather than to the ERGM.
EGP_fixedH<-function(family) family!=EGP_MS_DEGREES


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
      #Two independent reasons a cell can be unsupported, and they want
      #different notes: separability (the continuum STERGMs, under any
      #multi-toggle move) and a state-dependent move count (DS, under ~degrees
      #only -- see EGP_fixedH).
      nosep<-multi && !cap$multitog
      novH <-p=="DS" && !EGP_fixedH(fams[i])
      supported<-!nosep && !novH
      #DS is a special case for the engine too: its rate involves the size of
      #the move set, which for a multi-toggle family can only be had by
      #enumerating it, so thinning is unavailable there however bounded the
      #rate is.
      thin<-cap$thin && !(p=="DS" && multi)
      data.frame(process=p, constraint=cons[i],
                 supported=supported,
                 engine=if(!supported) "-" else if(thin) "thinning (exact)" else "enumeration",
                 note=if(nosep) "multi-toggle moves break formation/dissolution separability"
                      else if(novH) "|H| varies by state, so the equilibrium would be |H|exp(pot)"
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


#Validate the dyad-varying pacing factor.
#
#rate.esp multiplies each dyad's transition rate by exp(rate.esp * ESP), where
#ESP is the dyad's number of shared partners.  Because ESP of a dyad does not
#depend on the state of that dyad's own edge variable, the modulation is
#symmetric between the two states it connects, and so it leaves the ERGM
#equilibrium exactly where it was: it is a pure *opportunity* (search) effect,
#as distinct from changing a coefficient, which is a *preference* effect and
#does move the equilibrium.
#
#The symmetry argument is what constrains where this can be used, so the
#restrictions below are correctness requirements, not conveniences.
#
#rate.esp.cap caps the modulation at exp(rate.esp * min(ESP, cap)).  min(ESP,cap)
#is still a function of the dyad's neighbourhood alone, so the symmetry argument
#is untouched and the equilibrium is still exactly the specified ERGM.  What the
#cap buys is a tight thinning bound; see EGP_engine() and src/simEGP.c.
EGP_check_rate_esp<-function(rate.esp, rate.esp.cap=5, con, nw){
  if(is.null(rate.esp)||length(rate.esp)!=1L||!is.finite(rate.esp))
    stop("`rate.esp` must be a single finite number (0 disables it).")
  #Checked before the rate.esp==0 exit, so that a bad cap is never silently
  #accepted just because the modulation happens to be switched off.  Inf is
  #legal and means "no cap"; a *negative* cap is not merely odd but unsound,
  #since min(ESP,cap) would then be constant and negative, putting the true rate
  #above the bound the thinning engine assumes.
  if(is.null(rate.esp.cap)||length(rate.esp.cap)!=1L||!is.numeric(rate.esp.cap)||
     is.na(rate.esp.cap)||rate.esp.cap<0)
    stop("`rate.esp.cap` must be a single non-negative number, or Inf for no cap.\n",
         "  It is the largest shared-partner count the pacing factor responds to, so a\n",
         "  negative value has no meaning; it would also invalidate the rate bound the\n",
         "  thinning engine relies on, which assumes the modulation exponent lies in\n",
         "  [0, rate.esp.cap].")
  #A zero cap makes the modulation identically 1, so there is nothing left for
  #the restrictions below to protect.
  if(rate.esp==0||rate.esp.cap==0) return(invisible(NULL))
  if(is.directed(nw))
    stop("`rate.esp` is defined for undirected networks only.\n",
         "  Shared-partner counts have several inequivalent directed analogues (OTP, ITP,\n",
         "  OSP, ISP) and picking one silently would be a modelling decision, not a default.")
  if(con$ntoggles!=1L)
    stop("`rate.esp` cannot be combined with the constraint ",con$label,".\n",
         "  It is only equilibrium-preserving for single-dyad moves: a dyad's shared-partner\n",
         "  count is unchanged by toggling that dyad, which is what makes the rate\n",
         "  modulation symmetric.  A multi-dyad move (",con$ntoggles," toggles) alters the\n",
         "  neighbourhoods of the dyads involved, so ESP is no longer invariant, symmetry\n",
         "  fails, and the process would no longer have the specified ERGM equilibrium.\n",
         "  Use rate.esp with an unconstrained or dyad-level-constrained process.\n",
         "  (A symmetric generalisation to multi-dyad moves does exist -- the modulation\n",
         "  would have to be a function of the *unordered* set of toggled dyads and of the\n",
         "  rest of the graph, e.g. the mean shared-partner count of the dissolved and\n",
         "  formed dyads computed with both deleted -- but which one to use is a modelling\n",
         "  choice, so ergmgp does not pick one for you.)")
  invisible(NULL)
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
    #DS is refused where the legal-move count varies by state, because there its
    #equilibrium is not the ERGM the caller specified.  See EGP_fixedH above for
    #the counting argument and for why ~degrees is the only such case.
    if(proc=="DS" && !EGP_fixedH(spec$family))
      stop("Process DS cannot be simulated under ~",cntcon,".\n",
           "  DS's rate is A/(|H|exp(pot)), so its jump chain is uniform over the |H| legal\n",
           "  moves and its total exit rate is A exp(-pot).  Its equilibrium is therefore\n",
           "  proportional to |H|(x)exp(pot(x)), which is the requested ERGM only when |H| is\n",
           "  constant.  A tetrad is legal only if the two dyads it would form are absent, so\n",
           "  under ~",cntcon," |H| varies with the configuration and the equilibrium is tilted\n",
           "  towards states with more legal moves.  Simulating it would return draws from a\n",
           "  distribution that is not the one specified by `coef`.\n",
           "  DS is exact under every other constraint (|H| is fixed by the constraint\n",
           "  itself there), so use ~edges, ~odegrees, ~idegrees or a dyad-level constraint;\n",
           "  or use LERGM or CI, which are exact under ~",cntcon,".")
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
#  rate.esp, rate.esp.cap - the dyad-varying pacing factor, which bears on the
#            choice: the thinning bound has to dominate the modulation, and how
#            tightly it can do so depends on whether the cap is finite
#  verbose - logical; report the choice?
#
#Return value: 0L for enumeration, 1L for thinning.
EGP_engine<-function(con, proc, engine=c("auto","thinning","enumeration"),
                     rate.esp=0, rate.esp.cap=5, verbose=FALSE){
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
  #
  #A dyad-varying pacing factor changes that calculus, because the thinning
  #bound has to dominate the modulation over every dyad.  With a finite cap the
  #bound is exactly A*exp(rate.esp*rate.esp.cap): tight, state-independent, and
  #free of the per-proposal max-degree scan, so thinning is right even for the
  #unconstrained family (rate.esp is new, so there is no prior behaviour to
  #preserve there).  With no cap the bound must assume the most-connected dyad,
  #giving an acceptance probability of about exp(-rate.esp * max degree) --
  #unusable above a few hundred nodes -- so "auto" enumerates instead.
  esp.on<-rate.esp!=0 && rate.esp.cap>0
  use.thin<-switch(engine,
                   auto=if(esp.on) cap$thin && is.finite(rate.esp.cap)
                        else cap$thin && con$family!=EGP_MS_FREE,
                   thinning=TRUE,
                   enumeration=FALSE)
  if(use.thin && esp.on && !is.finite(rate.esp.cap))
    warning("engine=\"thinning\" with rate.esp=",rate.esp," and no cap: the rate bound has\n",
            "  to assume the most-connected dyad, so the acceptance probability is about\n",
            "  exp(-rate.esp * max degree).  This is still exact, but on any sizeable\n",
            "  network it may not finish.  Set `rate.esp.cap` (the default is 5).",
            call.=FALSE)
  if(use.thin && esp.on && rate.esp*rate.esp.cap>15)
    warning("rate.esp * rate.esp.cap = ",format(rate.esp*rate.esp.cap),
            " puts the thinning acceptance\n  probability below exp(-15), so the run may not",
            " finish; beyond about 700 the\n  inter-event time underflows to zero and the",
            " clock stops advancing altogether\n  (use.logtime=TRUE avoids that part).",
            "  Reduce `rate.esp.cap`.",call.=FALSE)
  if(verbose){
    message("EGP: process ",proc,", constraint ",con$label,", ",
            con$ntoggles," toggle(s) per event, engine ",
            if(use.thin) "thinning (exact uniformization)" else "enumeration",".")
    #The cap is a *model* parameter and its default is finite, so a run should
    #never be quietly capped: say what is in force and how to turn it off.
    if(esp.on)
      message("EGP: dyad-varying pacing, rate.esp=",rate.esp,
              if(is.finite(rate.esp.cap))
                paste0(", modulated by min(ESP, ",format(rate.esp.cap),
                       ").  The cap leaves the equilibrium unchanged but does change the ",
                       "dynamics; pass rate.esp.cap=Inf for the uncapped process.")
              else ", uncapped.")
    if(!use.thin && con$ntoggles>1L)
      message("EGP: note - enumerating a ",con$ntoggles,
              "-toggle move set is expensive; consider a process that supports thinning ",
              "(LERGM, CI, DS).  See EGPConstraintSupport().")
  }
  as.integer(use.thin)
}
