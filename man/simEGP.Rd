\name{simEGP}
\alias{simEGP}
\alias{simEGPTraj}

\title{
Simulate Trajectories from an ERGM Generating Process
}

\description{
Given an \code{\link[ergm]{ergm}} formula, simulate trajectories from a continuous time graph process having the specified ERGM as a limiting distribution.  A number of different processes are supported, and termination may be specified either by phenomenological time or event counts.
}

\usage{
simEGP(form, coef, events = 1, time = NULL, rate.factor = 1,
    time.offset = 0, event.offset = 0, process = c("LERGM", "CRSAOM",
    "CI", "DS", "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"),
    constraints = NULL, engine = c("auto", "thinning", "enumeration"),
    rate.esp = 0, rate.esp.cap = 5, use.logtime = FALSE, return.changetime = FALSE,
    changetime.offset = NULL, return.history = FALSE,
    return.networkDynamic = FALSE, verbose = TRUE,  trace.interval = 100,
    ...)

simEGPTraj(form, coef, events = 1, time = NULL, checkpoints = 1,
    rate.factor = 1, trajectories = 1, mc.cores = 1,
    log.sampling = FALSE, process = c("LERGM", "CRSAOM", "CI", "DS",
    "CDCSTERGM", "CFCSTERGM", "CSTERGM", "CTERGM"), constraints = NULL,
    engine = c("auto", "thinning", "enumeration"), rate.esp = 0,
    rate.esp.cap = 5, use.logtime = FALSE,
    return.changetime = FALSE, return.history = FALSE, verbose = TRUE,
    trace.interval = 100, statsonly = FALSE, monitor = NULL)
}

\arguments{
  \item{form}{
an \code{\link[ergm]{ergm}} formula defining terms for the EGP; the left-hand side must be a network object, whose properties are used to determine the state space.  For the \code{CSTERGM} process, a list containing two such formulas must be used, with named elements \code{formation} (for the formation model) and \code{dissolution} (for the dissolution model).
}
  \item{coef}{
vector of coefficients for the EGP; for the CSTERGMs, this should be a containing named elements \code{formation} and \code{dissolution}, each of which should be the coefficient vector for its respective model.
}
  \item{events}{
optionally, the number of simulated events to draw (if \code{time==NULL}); if \code{time} is specified, this is ignored.
}
  \item{time}{
optionally, the temporal length of the simulation; if not supplied, \code{events} is used instead to determine when to stop.
}
  \item{rate.factor}{
a multiplicative factor scaling the time evolution of the system; higher values correspond to faster dynamics.
}
  \item{time.offset}{
optionally, an initial \dQuote{clock} offset for the start of a trajectory; this allows time 0 (the start of the simulation interval) to be set to an arbitrary time point.  This is only used for book-keeping (e.g., when a trajectory is run as multiple segments), and does not affect e.g. the meaning of the \code{time} argument (which is always interpreted as units after the start time).
}
  \item{event.offset}{
optionally, an initial offset to the step or event count for the start of a trajectory (e.g., for trajectories being run in segments).  As with \code{time.offset}, this only affects book-keeping, and has no other effect.
}
  \item{process}{
the ERGM generating process to use (described below).
}
  \item{constraints}{
optionally, a one-sided formula restricting the sample space, using the same syntax as \code{\link[ergm]{ergm}}: any dyad-level constraint, or one of \code{~edges}, \code{~degrees}, \code{~odegrees}, \code{~idegrees}.  Not every combination of process and constraint is supported; see \sQuote{Constrained sample spaces} below, and \code{\link{EGPConstraintSupport}} for the full table.  Unsupported combinations raise an error rather than being silently ignored.
}
  \item{engine}{
how to simulate the process: \code{"thinning"} (exact uniformization; one change-statistic evaluation per candidate move), \code{"enumeration"} (walk the whole move set at every event), or \code{"auto"} (the default).  Both are exact; they differ only in cost.  Under \code{"auto"} the constrained families are thinned wherever the process admits a rate bound, while the unconstrained process stays on enumeration so that results predating the thinning engine are reproduced; and when \code{rate.esp} is non-zero, thinning is selected if and only if \code{rate.esp.cap} is finite, since only then is the bound tight enough to be useful.  See \sQuote{Constrained sample spaces} and \sQuote{Opportunity versus preference}.
}
  \item{rate.esp}{
a dyad-varying pacing factor: each edge variable's transition rate is multiplied by \code{exp(rate.esp * min(ESP, rate.esp.cap))}, where \code{ESP} is that dyad's number of shared partners.  Zero (the default) recovers the ordinary process.  Positive values make well-embedded dyads come up for reconsideration more often without changing how attractive they are, which leaves the equilibrium untouched; see \sQuote{Opportunity versus preference} below.  Undirected, single-dyad moves only.
}
  \item{rate.esp.cap}{
the largest shared-partner count the pacing factor responds to: dyads with more shared partners than this are paced as though they had exactly this many.  It leaves the equilibrium exactly where \code{rate.esp} alone leaves it, for the same reason, but it does change the dynamics, so it is a model parameter to be chosen and reported rather than a tuning constant to be adjusted until results settle.  It also makes exact simulation affordable: the transition rate is then bounded by \eqn{A e^{\lambda c}}, tightly, where an uncapped run has to assume the most-connected dyad and becomes unusable above a few hundred nodes.  Defaults to 5; pass \code{Inf} for the uncapped process.  Ignored when \code{rate.esp} is zero, though still validated.
}
  \item{use.logtime}{
logical; internally, use logarithmic timescale?  This can potentially protect against overflow or underflow when rates are extreme, but can reduce precision and adds some overhead.
}
  \item{return.changetime}{
logical; should we return a matrix with the last update times for each edge variable as a network attribute?
}
  \item{changetime.offset}{
optionally, an \code{n x n} matrix of last change times (for trajectories being resumed in process).
}
  \item{return.history}{
logical; return the entire event history as a network attribute?
}
  \item{return.networkDynamic}{
logical; retain the entire event history and return as a \code{\link[networkDynamic]{networkDynamic}} object?
}
  \item{verbose}{
logical; provide trace messages regarding simulation progress?
}
  \item{trace.interval}{
for verbose output, the interval at which messages should be printed (in events).
}
  \item{checkpoints}{
number of checkpoints at which the trajectory should be sampled (in addition to the initial state).
}
  \item{trajectories}{
number of independent trajectories to simulate (all start from the seed network, but evolve independently).
}
  \item{mc.cores}{
number of cores to use when simulating trajectories.
}
  \item{log.sampling}{
logical; should time points to sample be logarithmically spaced?
}
  \item{statsonly}{
logical; should only network statistics be retained (and not the graphs themselves)?
}
  \item{monitor}{
optionally, an \code{\link[ergm]{ergm}} formula with additional statistics to track.
}
  \item{\dots}{
additional arguments (currently unused).
}
}

\details{
An ERGM generating process (EGP) is a continuous time graph process with an equilibrium distribution having a known ERGM form.  See \code{\link{ergmgp}} for an overview of EGPs, including the specifications supported here.

\code{simEGP} generates a single trajectory from an EGP, with the EGP being specified via its graph potential (as a \code{\link[ergm]{ergm}} formula or pair thereof and associated coefficients) and its initial state being given by the left-hand side of the input formula.  The trajectory length can be specified either in terms of the number of transitions to be simulated (\code{events}) or the length of the trajectory in phenomenological time (\code{time}); only the latter leads to the specified ERGM equilibrium (since event times are not \dQuote{random} times, stopping after a fixed number of events biases the final state).  If desired for bookkeeping purposes, an offset can be added to the simulation clock (which otherwise starts at 0), event count (likewise), and most recent change times (also likewise).  By default, the return value is a \code{\link[network]{network}} object containing the final graph state, with network attributes giving the final time (\code{"Time"}), event count (\code{"Events"}), ERGM potential (\code{"Potential"}).  A square matrix containing the time of the most recent transition experienced by each edge variable can be returned as a network attribute (\code{"LastChangeTime"}) if \code{return.changetime} is selected.  By default, the entire event history is not stored (as it can become extremely large).  However, if \code{return.history==TRUE}, a matrix containing the event history is saved and returned as a network attribute (\code{"EventHistory"}). Alternately, setting \code{return.networkDynamic=TRUE} will lead to the event history being stored and the entire trajectory being returned as a \code{\link[networkDynamic]{networkDynamic}} object, with edge activities set based on the observed transitions.  This format may be easier to use for visualization, or to query the state of the network at an arbitrary point in the trajectory.  The \code{\link{durations}} function can be used to extract edge durations from such objects, as well.

For models with extreme transition rates, the option \code{use.logtime} may be useful for avoiding overflow or underflow; this only affects internal calculation, and not reported event times.  Note that logscale calculations add some overhead, and may be less precise in some cases than the default, so this option is not suggested unless specifically needed.  

\code{simEGPTraj} is a wrapper for \code{simEGP}, which adds additional capabilities for simulation of multiple trajectories and/or sampling of longer trajectories.  Each returned trajectory contains the initial state, as well as \code{simEGP} output from \code{checkpoints} points along the trajectory (including the end).  The default behavior (\code{checkpoints==1}) returns the initial and final states.  Checkpoints are evenly spaced (with termination criteria indicated as per \code{simEGP}) by default, or logarithmically spaced if \code{log.sampling==TRUE}.  Multiple independent trajectories can be simulated by setting \code{trajectories>1}; these can be run in parallel by setting \code{mc.cores>1}.  If desired, the model statistics can be returned without the graph state by choosing \code{statsonly=TRUE}, and a one-sided \code{monitor} formula can likewise be used to calculate additional statistics if desired (with similar functionality to the \code{ergm} \code{simulate} method).  Otherwise, \code{network.list} objects are returned containing the states in the respective trajectories.

Simulation itself follows the discrete event approach described in Butts (2023).  Transition hazards are computed for all edge variables (making the scaling no better then \eqn{O(N^2)} for each update, and are used to draw both the next event and the event time.  Because the cost of computing each transition is unrelated to waiting time, this algorithm can be quite efficient at simulating long time periods when events are sparse (unlike, e.g. a discrete-time algorithm that updates in every period).  By turns, however, trajectories can become quite expensive (per unit phenomenological time) when event rates are high.  This issue is especially pronounced for the \code{CRSAOM} and \code{DS} processes, which can both generate very high transition rates in some cases.  Unless otherwise specified, event histories are not stored, so storage costs are by default unrelated to trajectory length.  Care should be taken when storing event histories, as they can become quite large when transition rates are high. 

To obtain equilibrium graph distributions from an EGP, it is generally (much) more efficient to use the \code{simulate} functions in the \code{\link[ergm]{ergm}} package: they employ MCMC algorithms that are unconstrained by the need to follow realistic trajectories, and that are optimized for rapid mixing.  (In particular, note that many systems can become \emph{kinetically trapped}, spending very long periods in metastable states that are far from equilibrium.  This can be a real-world phenomenon, but is not always desirable from a computational point of view.  Functions such as \code{simEGP} are intended to faithfully reproduce such dynamics, while MCMC algorithms are intended to avoid them.)  Comparison of late-phase draws from a \code{simERGMPot} trajectory with equilibrium ERGM draws can be used to evaluate convergence to equilibrium behavior (where desired); alternately, \code{simEGP} can be seeded with ERGM draws to follow trajectories from equilibrated states.  Consult the \code{ergm} package documentation for details.
}

\section{Constrained sample spaces}{
Supplying \code{constraints} restricts the process to a subset of the graph space.  Because an EGP moves in continuous time, this is not merely a matter of rejecting illegal states: it changes what an \emph{event} is.

The unconstrained processes of Butts (2023) move by toggling a single dyad.  Under a dyad-level constraint (\code{blocks}, \code{observed}, \code{fixedas}, \code{fixallbut}, \code{egocentric}, and combinations) that is still true -- the fixed dyads are simply removed from consideration -- so those constraints are supported for every process.

The \emph{count-preserving} constraints are different: under \code{~edges} or a degree constraint no single toggle is legal at all, since every toggle changes the edge count and both its endpoints' degrees.  A move therefore has to become a \emph{set} of simultaneous toggles.  ergmgp implements three shapes:
\describe{
  \item{\code{~edges}}{two toggles: dissolve one edge, form one non-edge, anywhere in the graph.  The edge count is held fixed and the degree distribution is free to move.  Available for directed and undirected networks.  This is the same move as \code{ergm}'s \code{ConstantEdges} proposal.}
  \item{\code{~odegrees}, \code{~idegrees}}{two toggles: rewire an edge while holding one endpoint fixed.  Directed networks only.}
  \item{\code{~degrees}}{four toggles: the tetrad or double-edge swap, mirroring \code{ergm}'s \code{CondDegreeTetrad} proposal.  Undirected networks only.}
}
At most one count-preserving constraint can be in force, since each dictates a different move shape; combining two is an error.  (\code{~edges + degrees} is not such a case: \code{~degrees} fixes every degree and hence their sum, and \code{ergm} drops the redundant term before ergmgp sees it.)

The resulting process has the constrained ERGM as its equilibrium provided the rate function satisfies the usual detailed-balance form applied to the potential difference of the whole move, and the move graph is symmetric and connected on the constrained space.  This is an extension of the model class, not just of the implementation.

Which of these to use is a modelling decision and not an innocuous one.  \code{~degrees} says each node has a fixed tie budget, so only the \emph{composition} of ties can change; the unconstrained process lets volume and composition move together, so a rise in any composition statistic is confounded with a rise in the number of ties; \code{~edges} sits between them, fixing total volume while allowing ties to redistribute across nodes -- the constraint to use when the hypothesis is about \emph{concentration} rather than composition.

The continuum STERGMs (\code{CSTERGM}, \code{CDCSTERGM}, \code{CFCSTERGM}) cannot be extended this way: any of these moves is simultaneously a formation and a dissolution, so the potential difference admits no principled split into formation and dissolution parts and separability breaks down.  Requesting a count-preserving constraint for these processes is an error.

Two exact engines are available, differing only in cost:
\describe{
  \item{\code{"enumeration"}}{walks the whole legal move set at each event, exactly as the unconstrained algorithm walks all dyads.  Works for every process, but the move set has \eqn{O(E^2)} members under \code{~degrees} and \eqn{E(D-E)} under \code{~edges}, and multi-toggle change statistics must be obtained by applying and rolling back the toggles (terms such as \code{gwesp} supply no joint change statistic), so this is one to two orders of magnitude slower than the unconstrained case.}
  \item{\code{"thinning"}}{uniformization.  Where the transition rate is bounded above by a constant \eqn{\Lambda} -- true for \code{LERGM} and \code{CI} (\eqn{\Lambda = A}), and trivially for \code{DS}, whose rate does not depend on the move -- candidates are drawn uniformly from a fixed superset of the move set, and accepted with probability \eqn{rate/\Lambda}.  Rejections are self-loops: time advances, the state does not.  This is \emph{exact}, not an approximation, and costs one change-statistic evaluation per candidate instead of one per member of the move set.  It is the default wherever it applies.}
}
\code{DS} under a count-preserving constraint is restricted to enumeration, because its rate involves the size of the move set, which can only be obtained by enumerating it; enumeration is cheap for \code{DS} since no change statistics are needed during the sweep.

That same dependence on the move count makes \code{DS} the one process for which a count-preserving constraint is not merely a matter of cost, and \code{DS} under \code{~degrees} is refused.  Its rate is \eqn{A/(|H|\exp(pot))}, so the total exit rate is \eqn{A\exp(-pot)} and the embedded jump chain is uniform over the \eqn{|H|} legal moves.  The move graph is symmetric, so that jump chain's stationary distribution is proportional to \eqn{|H|(x)}, and \eqn{\pi = \nu/q} for a continuous-time chain; the equilibrium is therefore proportional to \eqn{|H|(x)\exp(pot(x))}, which is the requested ERGM only where \eqn{|H|} is constant.  It is constant for the single-toggle families, for \code{~edges} (\eqn{E(D-E)} moves) and for the rewire families (\eqn{\sum_t d_t(n-1-d_t)}), each pinned by what the constraint itself pins, and it stays constant when these are composed with a dyad-level constraint.  Under \code{~degrees} it is not: a tetrad is legal only when the two dyads it would form are currently absent, so the count depends on the configuration rather than on the degree sequence.  \code{DS} would then equilibrate to a distribution tilted towards states offering more legal moves, and the combination raises an error rather than returning draws from a model other than the one requested.  \code{LERGM} and \code{CI} are unaffected, their rates depending on the move and not on how many moves there are.

The size of that tilt is computable rather than merely bounded, since \eqn{\pi_{DS}} is exactly \eqn{\pi} reweighted by \eqn{|H|}: the induced shift in a statistic \eqn{g} is \eqn{E[g|H|]/E[|H|] - E[g]} under the target.  It is small and falls quickly with size -- at fixed mean degree, roughly \eqn{n^{-1.1}} absolutely and \eqn{n^{-1.65}} relative to the statistic, about 0.09\% at \eqn{n = 14} and 0.001\% at \eqn{n = 200}, and smaller again at higher density.  The refusal is therefore not a claim that the error is large; it is that the process does not have the stated equilibrium, and that the discrepancy is worst in precisely the small, sparse regime where an equilibrium is most likely to be checked exactly.

The realised constraint and engine are recorded on the returned network as the \code{"Constraint"} and \code{"Engine"} attributes; for the thinning engine, \code{"Acceptance"} gives the fraction of candidates accepted.  A low acceptance rate costs speed but never correctness.  A modulated run additionally records \code{"RateESP"} and \code{"RateESPCap"}, so that a network simulated with \code{verbose=FALSE} still carries the pacing model it was produced under rather than leaving a defaulted cap implicit.

Use \code{EGPConstraintSupport} to print the full process-by-constraint table.  Any unsupported combination is an error: this function no longer accepts and ignores unrecognised arguments.
}

\section{Opportunity versus preference}{
\code{rate.factor} is a scalar pacing constant: it sets the timescale, but is the same for every edge variable.  \code{rate.esp} generalises it to a \emph{dyad-varying} pacing factor,
\deqn{A_{th} = A \exp(\lambda\, \mathrm{ESP}_{th}),}
where \eqn{\mathrm{ESP}_{th}} is the number of partners \eqn{t} and \eqn{h} share and \eqn{\lambda} is \code{rate.esp}.  Well-embedded dyads are then \emph{reconsidered} more often, without any change to how attractive a tie between them is once considered.

\strong{This does not move the equilibrium.}  A modulation \eqn{s} of the rates preserves detailed balance whenever \eqn{s(x,y) = s(y,x)}, and \eqn{\mathrm{ESP}_{th}} is unchanged by toggling \eqn{(t,h)} -- a dyad is not part of its own shared-partner count -- so the modulation applies identically in both directions:
\deqn{\pi(x) r'(x \to y) = s\,\pi(x) r(x \to y) = s\,\pi(y) r(y \to x) = \pi(y) r'(y \to x).}
The limiting distribution is therefore exactly the ERGM implied by \code{form} and \code{coef}, for any \eqn{\lambda}.

That gives a clean separation between two things an ERGM cannot tell apart:
\describe{
  \item{\emph{opportunity}}{a symmetric rate modulation (\code{rate.esp}).  Changes which ties are easy to \emph{find}.  The equilibrium is invariant; only the dynamics and the timescale move.}
  \item{\emph{preference}}{a change to \code{coef}.  Changes how attractive a tie is.  The equilibrium moves.}
}
The motivating case is friend-recommendation systems (\dQuote{people you may know}), which surface friends-of-friends: they alter the choice set rather than the utility.  Modelled as \code{rate.esp}, such a system provably cannot change the equilibrium level of any network statistic, however large \eqn{\lambda} is -- so any effect it has is a transient one, invisible to equilibrium analysis.  Modelled as an increase in a closure coefficient, it changes the equilibrium instead.  Which of the two is the better description is an empirical question, and the two make different predictions.

Restrictions, both of which are correctness requirements rather than conveniences.  \code{rate.esp} is available only for undirected networks, because shared-partner counts have several inequivalent directed analogues and choosing one silently would be a modelling decision.  And it is available only for single-dyad moves: under a count-preserving constraint a move toggles several dyads at once, altering the neighbourhoods of the dyads involved, so \eqn{\mathrm{ESP}} is no longer invariant, the symmetry argument fails, and the equilibrium would silently cease to be the specified ERGM.  Both cases raise an error.

Note that this rules out \code{rate.esp} in combination with \code{~edges} as well as with \code{~degrees}: the edge dissolved by an \code{~edges} swap may itself be part of the formed dyad's shared-partner count, so the modulation differs between the two states the move connects.  A symmetric generalisation does exist -- the modulation would have to be a function of the \emph{unordered} pair of toggled dyads and of the rest of the graph, for instance the mean shared-partner count of the dissolved and formed dyads computed with both deleted -- but there are several such functions and they are not equivalent, so ergmgp does not choose one on the user's behalf.  A study wanting both a fixed edge count and an opportunity effect must therefore either specify that choice itself or use a dyad-level constraint instead.

\subsection{The cap, and why it is not optional at scale}{
Exact simulation by thinning needs an upper bound on the transition rate, and the modulation \eqn{\exp(\lambda\,\mathrm{ESP}_{th})} has no useful bound: all that can be said in general is \eqn{\mathrm{ESP}_{th} \le \max_v \deg(v)}.  On a real network that is fatal rather than merely loose -- with a maximum degree in the hundreds and \eqn{\lambda = 0.5} the acceptance probability is around \eqn{e^{-100}}, which is not slow but zero -- so an uncapped run must enumerate, and enumeration does not scale either.  Uncapped, the feature is unusable much above \eqn{n = 200}.

\code{rate.esp.cap} caps the modulation at
\deqn{A_{th} = A \exp(\lambda \min(\mathrm{ESP}_{th}, c)),}
which bounds the rate by exactly \eqn{A e^{\lambda c}}.  That bound is tight, does not depend on the state, and costs nothing to evaluate, so thinning becomes exact \emph{and} affordable: \code{engine="auto"} selects thinning when the cap is finite and enumeration when it is not.  As a rule of thumb the acceptance probability scales like \eqn{e^{-\lambda c}}, so \eqn{\lambda c \lesssim 3} is comfortable, \eqn{\lambda c \approx 5}--\eqn{10} is expensive, and beyond \eqn{\lambda c \approx 15} a run may not finish; far beyond it the inter-event time underflows to zero and the clock stops advancing altogether, which \code{use.logtime=TRUE} avoids.  Both of the outer cases warn.

Capping does not disturb any of the above.  \eqn{\min(\mathrm{ESP}_{th}, c)} is still a function of the dyad's neighbourhood excluding the dyad itself, so it is still unchanged by toggling \eqn{(t,h)}, the modulation is still symmetric, and the equilibrium is still exactly the ERGM implied by \code{coef}.  What the cap does change is the dynamics, which is why it is a \emph{model parameter and not a tuning constant}: it should be chosen and reported, not raised until results stop moving.  There is a substantive case for it as the better model, too -- a recommendation surface has finitely many slots, and the boost from a sixth mutual friend is not the boost from a second, so saturation is the more defensible functional form and the uncapped exponential was never especially plausible at high \eqn{\mathrm{ESP}}.

The default is \code{rate.esp.cap = 5}.  Pass \code{Inf} for the uncapped process.
}

One further practical point: \code{rate.esp} rescales the clock.  Raising every embedded dyad's rate raises the total event rate by roughly \eqn{\exp(\lambda \bar{L})}, where \eqn{\bar{L}} is the mean shared-partner count, so a run to a given \code{time} draws proportionately more events and costs proportionately more.  The cap bounds that too, at \eqn{e^{\lambda c}}, which is the other reason a capped run is affordable at a horizon an uncapped one is not.  Comparing a modulated run against an unmodulated one at equal \code{time} is statistically the right comparison -- both approach the same equilibrium, the modulated one faster -- but the modulated arm may do an order of magnitude more work.  Budget accordingly, or compare at equal \code{events} when only the embedded jump chain is of interest.
}

\value{
  For \code{simEGP}, a network object containing the final graph state, with network attributes \code{Time}, \code{Events}, and \code{Potential} listing the time, event count, and ERGM potential at the end of the simulation interval.  See above for additional attributes that may be added if history retention is activated.  If \code{return.networkDynamic==TRUE}, then the return value is instead a \code{\link[networkDynamic]{networkDynamic}} object containing the event history as edge activity data; be aware that an edge will exist in this object if any corresponding edge is ever active, so the raw graph state should not be used to access the final system state.  Instead, use the \code{\link[networkDynamic]{network.extract}} method to query the network state at the desired time point.
  
  For \code{simERGMPotTraj}, a list containing the simulated trajectories.  These are either \emph{network.list} objects containing the networks at each checkpoint (with time, step, and potential attribute as described above), or else matrices of trace statistics (if \code{statsonly==TRUE}).  Note that the statistics are in any event included as an attribute to each network list, so the effect of \code{statsonly==TRUE} is simply not to retain the graph states.
}

\references{
Butts, Carter T.  (2023).  \dQuote{Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.} \emph{Journal of Mathematical Sociology}.  \doi{10.1080/0022250X.2023.2180001}
}

\author{
Carter T. Butts \email{buttsc@uci.edu}
}

\note{
Using \code{steps} to control trajectory termination will lead to biased samples (sometimes severely so); this is because transitions are not random times. If your goal is to obtain equilibrium draws (or draws en route thereto), use \code{time} to set the stopping point.  See \code{\link{EGPRateEst}} for a simple tool for calibrating simulation times.
}

\seealso{
\code{\link{ergmgp}} for information on EGPs, \code{\link[ergm]{ergm}} for information on ERGM specifications, \code{\link{EGPHazard}}, \code{\link{EGPRateEst}}, \code{\link[networkDynamic]{networkDynamic}}
}

\examples{
\donttest{
#Small example of 2-ribbon generation
n<-100
set.seed(1331)
net<-network.initialize(n,directed=FALSE)
sim<-simEGP(net~edges+kstar(2)+nsp(1:2), 
    coef=c(109-log(n),-25,-1.25,3.25), time=100, process="LERGM",
    verbose = TRUE)
plot(sim) #Return value is a single network

#Generate a trajectory showing the process at multiple stages
set.seed(1331)
sim<-simEGPTraj(net~edges+kstar(2)+nsp(1:2), 
    coef=c(109-log(n),-25,-1.25,3.25), time=100, checkpoints = 5, 
    trajectories = 2, mc.cores = 1, log.sampling = TRUE, 
    process = "LERGM", verbose = TRUE) 
length(sim)==2        #One entry per simulated trajectory
op<-par(mfrow=c(2,3))
for(i in 1:6)         #Show the first trajectory
  plot(sim[[1]][[i]],main=paste("Time",round(sim[[1]][[i]]\%n\%"Time",2)))
summary(sim[[2]]~edges+kstar(2))  #Show selected stats from the second
attributes(sim[[1]])  #Show precomputed statistics
par(op)

#A simple example with statsonly
set.seed(1331)
sim<-simEGPTraj(net~edges+esp(0), coef = c(log(2)-log(n), -1), time = 200,
    checkpoints = 25, process = "LERGM", statsonly = TRUE, 
    monitor = ~triangle)
sim                   #Note the monitor stat
op<-par(mfrow=c(1,1))
plot(sim[,"Time"], sim[,"edges"], type = "l")  #Time by edge count
lines(sim[,"Time"], sim[,"esp0"], col = 2)     #Add ESP(0)s
par(op)
}
}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ graphs }% __ONLY ONE__ keyword per line

