\name{EGPConstraintSupport}
\alias{EGPConstraintSupport}

\title{
Which Constraints Each ERGM Generating Process Supports
}

\description{
Print the table of process/constraint combinations that \code{\link{simEGP}} can simulate, which engine each uses, and why the unsupported combinations are unsupported.
}

\usage{
EGPConstraintSupport(process = NULL)
}

\arguments{
  \item{process}{
optionally, the name of a single process to report on; by default all eight are shown.
}
}

\details{
Constraining an EGP is not simply a matter of rejecting illegal states, because it changes what an \emph{event} is.  The unconstrained processes of Butts (2023) move by toggling a single dyad.  A dyad-level constraint leaves that intact, merely removing the fixed dyads from consideration, and is therefore available for every process.  A count-preserving constraint does not: no single toggle preserves the edge count, let alone a degree sequence, so a move must become a set of simultaneous toggles -- two for \code{~edges} (dissolve one edge, form one non-edge) and for \code{~odegrees} and \code{~idegrees} (rewire an edge about a fixed endpoint), four for \code{~degrees} (the tetrad).  At most one such constraint can be in force at a time.

Three consequences show up in this table.  First, the continuum STERGMs (\code{CSTERGM}, \code{CDCSTERGM}, \code{CFCSTERGM}) cannot accept count-preserving constraints at all: every one of these moves is simultaneously a formation and a dissolution, so the potential difference admits no principled split into formation and dissolution parts.  Second, the choice of engine depends on whether the process has a bounded transition rate.  Where it does -- \code{LERGM} and \code{CI}, whose rates are bounded by the pacing factor \eqn{A} -- the exact thinning (uniformization) engine can be used, which costs one change-statistic evaluation per candidate move instead of one per member of the move set.  Where it does not, the move set must be enumerated at every event, which is exact but slow.  \code{DS} is a special case: its rate is bounded, but the bound involves the size of the move set, which for a multi-toggle family can only be had by enumerating it -- so \code{DS} enumerates too, though cheaply, since a move-independent rate needs no change statistics during the sweep.

Third, \code{DS} is refused under \code{~degrees}, and this one is a matter of correctness rather than cost.  The differential-stability rate is \eqn{A/(|H|\exp(pot))}, where \eqn{|H|} is the number of legal moves, so the total exit rate is \eqn{A\exp(-pot)} and the embedded jump chain is uniform over those moves.  The move graph is symmetric, so that jump chain has stationary distribution proportional to \eqn{|H|(x)}, and a continuous-time chain has \eqn{\pi = \nu/q}; hence \code{DS} equilibrates to \eqn{|H|(x)\exp(pot(x))}, which is the requested ERGM exactly when \eqn{|H|} is constant.  It is constant under every constraint here but one: the move counts are \eqn{E(D-E)} for \code{~edges} and \eqn{\sum_t d_t(n-1-d_t)} for the rewire families, each fixed by what the constraint itself fixes.  A tetrad, however, is legal only if the two dyads it would form are currently absent -- a fact about the configuration, not the degree sequence -- so under \code{~degrees} the count varies and the equilibrium is tilted towards states offering more moves.  Rather than return draws from a distribution other than the one \code{coef} specifies, \code{\link{simEGP}} refuses the combination.  Use \code{LERGM} or \code{CI} under \code{~degrees}, or \code{DS} under one of the other constraints.

Dyad-level constraints are recognised structurally, as those exposing a free-dyad set: \code{blocks}, \code{observed}, \code{fixedas}, \code{fixallbut}, \code{egocentric} and combinations thereof, together with any future \code{\link[ergm]{ergm}} constraint of the same kind.  They compose with the count-preserving constraints by intersection.
}

\value{
Invisibly, a data frame with one row per process/constraint pair and columns \code{process}, \code{constraint}, \code{supported}, \code{engine}, and \code{note}.  The table is also printed.
}

\references{
Butts, C.T.  (2023).  \dQuote{Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.}  \emph{Journal of Mathematical Sociology}.
}

\author{
Carter T. Butts \email{buttsc@uci.edu}
}

\seealso{
\code{\link{simEGP}}
}

\examples{
#What can a longitudinal ERGM be constrained by?
EGPConstraintSupport("LERGM")

#...and why can a continuum STERGM not be constrained by degree?
EGPConstraintSupport("CSTERGM")
}

\keyword{graphs}
