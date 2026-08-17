# ergmgp 0.1-2.9000

## Constrained sample spaces

* `simEGP()` and `simEGPTraj()` gain a `constraints` argument, using `ergm`'s
  formula syntax. Supported: dyad-level constraints (`blocks`, `observed`,
  `fixedas`, `fixallbut`, `egocentric`, and any other constraint exposing a
  free-dyad set) for every process, and the count-preserving constraints
  `~edges` (either directedness), `~degrees` (undirected), `~odegrees` and
  `~idegrees` (directed) for every process except the continuum STERGMs.

  Count-preserving constraints require an event to toggle several dyads at
  once (two for an `~edges` swap or a rewire, four for a tetrad), which extends
  the process class of Butts (2023) rather than just its implementation; see the
  *Constrained sample spaces* section of `?simEGP`.

* `~edges` holds the total number of ties fixed while leaving the degree
  distribution free: a move dissolves one edge and forms one non-edge anywhere
  in the graph, the same move as `ergm`'s `ConstantEdges` proposal. It fills the
  gap between `~degrees`, which fixes every node's tie budget and so isolates
  the *composition* of ties, and the unconstrained process, in which volume and
  composition move together and a rise in any composition statistic is
  confounded with a rise in the number of ties. `~edges` is the constraint to
  use when the hypothesis concerns the *concentration* of a fixed volume of ties
  — whether ties pile up on already-popular nodes — which neither of the others
  can express.

  Only one count-preserving constraint may be in force, since each dictates a
  different move shape. `~edges + degrees` is not such a conflict: fixing every
  degree fixes their sum, and `ergm` drops the redundant term before ergmgp sees
  it.

* New `EGPConstraintSupport()` prints the full process-by-constraint table,
  saying for each combination whether it is supported, which engine is used,
  and why anything unsupported is unsupported.

* New `engine` argument selects between two exact simulation engines:
  `"enumeration"` (walk the whole move set at each event; works for every
  process) and `"thinning"` (uniformization; one change-statistic evaluation
  per candidate move). Thinning is exact, not approximate, and applies where
  the transition rate is bounded — `LERGM`, `CI`, and `DS`. On a typical
  workload (`edges + nodematch + gwesp`, n = 50, LERGM) simulating under
  `~degrees` runs at roughly 12,000 events/sec thinned versus 4 events/sec
  enumerated.

  `engine = "auto"` (the default) uses thinning for constrained processes and
  leaves the *unconstrained* process on enumeration, so existing results are
  unaffected. Unconstrained simulation can opt in with `engine = "thinning"`,
  which is roughly 28× faster on the same workload.

* Simulated networks now carry `"Constraint"` and `"Engine"` attributes, plus
  `"Acceptance"` when thinning was used.

## Behaviour changes

* **`simEGP()` no longer silently ignores unrecognised arguments.** Its `...`
  was previously documented as "not currently used", which meant any
  misspelled or unsupported option was accepted and quietly dropped — a call
  passing `constraint = ~odegrees` ran happily *unconstrained* and returned
  results that looked plausible but answered a different question. Unknown
  arguments are now an error. (`constraint`, singular, partial-matches the new
  `constraints` argument and behaves as intended.)

* Unsupported process/constraint combinations raise an error naming the
  process, the constraint, and the reason, rather than being ignored. In
  particular, the continuum STERGMs (`CSTERGM`, `CDCSTERGM`, `CFCSTERGM`)
  refuse count-preserving constraints: every such move is simultaneously a
  formation and a dissolution, so the potential difference admits no principled
  split into formation and dissolution parts.

* A constraint that leaves no free dyads is an error rather than a simulation
  in which no event can ever occur.

## Internals

* New `src/moveset.{h,c}` defines the move sets and exposes each both as an
  enumeration cursor and as a uniform sampler over a fixed superset. The rate
  algebra in `src/simEGP.c` is unchanged, factored into `score_move()`.

* Multi-dyad change statistics go through `ergm`'s `ChangeStats()`, which is
  required rather than merely convenient: terms such as `gwesp` supply only a
  `c_function`, so a joint change statistic can only be obtained by applying
  and rolling back the toggles.

* First test suite for the package (`tests/testthat/`), covering constraint
  invariance, dyad-level restriction, the explicit-error behaviour, agreement
  between the two engines on mean inter-event time, and agreement of the
  constrained equilibrium with `ergm::simulate()` under both `~degrees` and
  `~edges`. The `~edges` arm additionally monitors a degree count, which is
  exactly the marginal `~degrees` pins and `~edges` frees, and so the one a
  wrong 2-toggle move set could get wrong while still reproducing the model
  terms.
