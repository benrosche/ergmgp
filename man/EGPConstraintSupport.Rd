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
Constraining an EGP is not simply a matter of rejecting illegal states, because it changes what an \emph{event} is.  The unconstrained processes of Butts (2023) move by toggling a single dyad.  A dyad-level constraint leaves that intact, merely removing the fixed dyads from consideration, and is therefore available for every process.  A degree-preserving constraint does not: no single toggle preserves a degree sequence, so a move must become a set of simultaneous toggles -- two for \code{~odegrees} and \code{~idegrees}, four for \code{~degrees}.

Two consequences show up in this table.  First, the continuum STERGMs (\code{CSTERGM}, \code{CDCSTERGM}, \code{CFCSTERGM}) cannot accept degree constraints at all: a rewire is simultaneously a formation and a dissolution, so the potential difference admits no principled split into formation and dissolution parts.  Second, the choice of engine depends on whether the process has a bounded transition rate.  Where it does -- \code{LERGM} and \code{CI}, whose rates are bounded by the pacing factor \eqn{A}, and \code{DS}, whose rate does not depend on the move at all -- the exact thinning (uniformization) engine can be used, which costs one change-statistic evaluation per candidate move instead of one per member of the move set.  Where it does not, the move set must be enumerated at every event, which is exact but slow.

Dyad-level constraints are recognised structurally, as those exposing a free-dyad set: \code{blocks}, \code{observed}, \code{fixedas}, \code{fixallbut}, \code{egocentric} and combinations thereof, together with any future \code{\link[ergm]{ergm}} constraint of the same kind.  They compose with the degree constraints by intersection.
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
