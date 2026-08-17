# ergmgp 0.1-2.9000

## Dyad-varying pacing: separating opportunity from preference

* `simEGP()` and `simEGPTraj()` gain `rate.esp`, which generalises the scalar
  pacing constant `rate.factor` to a dyad-varying one:

      A_th = rate.factor * exp(rate.esp * ESP_th)

  where `ESP_th` is the number of partners the dyad shares. Dyads with more
  mutual contacts are *reconsidered* more often, without any change to how
  attractive a tie between them is once considered.

* **This provably does not move the equilibrium.** A modulation of the rates
  preserves detailed balance whenever it is symmetric in the two states a
  toggle connects, and a dyad's shared-partner count is unchanged by toggling
  that dyad, so the modulation applies identically in both directions. The
  limiting distribution remains exactly the ERGM implied by `form` and `coef`,
  for any value of `rate.esp`.

  This separates two mechanisms an ERGM cannot tell apart: a symmetric rate
  modulation is an *opportunity* effect (which ties are easy to find; equilibrium
  invariant, dynamics change), whereas a change to `coef` is a *preference*
  effect (how attractive a tie is; equilibrium moves). The motivating case is
  friend-recommendation systems, which alter the choice set rather than the
  utility, and which therefore -- on this reading -- can have no equilibrium
  effect at all, only a transient one.

* `rate.esp` is restricted to undirected networks (directed shared-partner
  counts have several inequivalent definitions) and to single-dyad moves (under
  a count-preserving constraint a move toggles several dyads, ESP is no longer
  invariant, and the symmetry argument fails). Both cases are errors rather than
  silent approximations.

  The second restriction covers `~edges` as well as the degree constraints: the
  edge an `~edges` swap dissolves may itself be part of the formed dyad's
  shared-partner count, so the modulation differs between the two states the
  move connects. A symmetric generalisation exists — the modulation would have
  to be a function of the *unordered* pair of toggled dyads and of the rest of
  the graph, e.g. the mean shared-partner count of the dissolved and formed
  dyads computed with both deleted — but several inequivalent such functions
  exist, so ergmgp does not choose one for you. A design wanting both a fixed
  edge count and an opportunity effect has to specify that choice itself.

* **New `rate.esp.cap`, defaulting to 5, caps the modulation** at
  `exp(rate.esp * min(ESP, cap))`. Without it the feature does not scale, and
  not in the ordinary sense: exact simulation by thinning needs an upper bound
  on the rate, the only general bound on ESP is the maximum degree, and on a
  real network with `rate.esp = 0.5` that puts the acceptance probability
  around `e^-100`. Not slow — zero. So an uncapped run has to enumerate, and
  enumeration does not scale either; the feature was unusable much above
  n = 200.

  With a finite cap the bound is exactly `A*exp(rate.esp*cap)` — tight,
  state-independent, and free of the per-proposal max-degree scan — so
  `engine="auto"` now selects thinning when the cap is finite and enumeration
  when it is not. Acceptance scales roughly as `exp(-rate.esp*cap)`; ergmgp
  warns when that product exceeds 15, and when thinning is requested with no
  cap at all.

  **The invariance result is untouched.** `min(ESP, cap)` is still a function of
  the dyad's neighbourhood excluding the dyad itself, so it is still unchanged
  by toggling that dyad, the modulation is still symmetric, and the equilibrium
  is still exactly the ERGM implied by `coef`. (The tempting refinement
  `min(ESP, min(cap, max degree))` is *not* safe, and the code says so: maximum
  degree changes when the dyad is toggled, so that version would break symmetry
  silently, while still producing plausible-looking networks.)

* **The cap's default is finite, so it changes results** for an existing
  `rate.esp` call. The equilibrium is invariant either way, but the transient is
  not — and the transient is the entire content of an opportunity effect. It is
  therefore announced rather than assumed: a run reports the cap in force under
  `verbose`, and every modulated network records `"RateESP"` and `"RateESPCap"`
  so that a `verbose=FALSE` run still carries the model it was produced under.
  `rate.esp.cap = Inf` restores the uncapped process exactly.

  The cap is a model parameter, not a numerical tolerance. There is a
  substantive case for it: a recommendation surface has finitely many slots, so
  the boost from a sixth mutual friend is not the boost from a second, and
  saturation is a more defensible functional form than an unbounded exponential.

* Fixed: the `engine="auto"` fallback for a non-zero `rate.esp` tested
  `identical(engine, c("auto","thinning","enumeration"))`, i.e. the *unevaluated
  default*, so it never fired for an explicit `engine="auto"` nor for anything
  arriving from `simEGPTraj()`. Under a dyad-level constraint, where `auto`
  gives thinning, those paths silently used the loose bound. The rule now lives
  in `EGP_engine()` and applies however `auto` is reached.

* Fixed: `simEGPTraj()` did not validate `rate.esp` up front, so a misuse
  surfaced from inside a worker — and on Unix `mclapply()` returns `try-error`
  objects rather than stopping, so it could surface as a malformed return value
  rather than as an error.

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

  `rate.esp` remains unavailable under `~edges`, for the same reason it is
  unavailable under `~degrees` — see below.

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

* **`DS` refuses `~degrees`**, and unlike the other refusals this one is about
  correctness rather than about what the code can express. The
  differential-stability rate is `A/(|H| exp(pot))`, so the total exit rate is
  `A exp(-pot)` and the embedded jump chain is uniform over the `|H|` legal
  moves. The move graph is symmetric, so that jump chain's stationary
  distribution is proportional to `|H|(x)`, and `pi = nu/q` for a
  continuous-time chain — hence `DS` equilibrates to `|H|(x) exp(pot(x))`,
  which is the requested ERGM exactly when `|H|` is constant.

  It is constant almost everywhere: `E(D-E)` moves under `~edges`,
  `sum_t d_t(n-1-d_t)` under the rewire constraints, the dyad count or the free
  dyad count for the single-toggle families — each pinned by whatever the
  constraint itself pins, and unchanged by composing with a dyad-level
  constraint, since frozen dyads never move. `~degrees` is the exception, and
  for a reason specific to the tetrad: it is legal only when the two dyads it
  would *form* are currently absent, which is a fact about the configuration
  and not about the degree sequence. On a 14-node, 28-edge model the count
  ranges over 250–310 across equilibrium draws, tilting the equilibrium towards
  states that offer more moves. Rather than return draws from a model other
  than the one `coef` specifies, the combination is an error. `LERGM` and `CI`
  are unaffected: their rates depend on the move, not on how many moves exist.

  The tilt is small, and shrinks fast: since it is exactly a reweighting by
  `|H|`, its size can be computed rather than simulated, and at fixed mean
  degree the shift in a gwesp statistic falls as roughly `n^-1.1` in absolute
  terms and `n^-1.65` relative — 0.09% of the statistic at n=14, 0.001% at
  n=200 — with raising the density damping it further. `|H|` is a count of
  order `E^2` whose random part is only the blocked re-pairings, of order `E^2`
  times the density, so its relative fluctuation dies quickly. This is a reason
  to state the magnitude, not a reason to allow the combination: the process
  does not have the equilibrium the package says it has, and the discrepancy is
  largest in exactly the small, sparse regime where a user is most likely to be
  checking an equilibrium against an exact calculation.

  There is no version of this that keeps `DS` as `DS`. Dropping the `1/|H|`
  restores detailed balance but destroys the property that defines the
  process — the exit rate would then depend on the number of available moves
  rather than on the potential alone.

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

  A separate test covers `DS` under `~edges`, the constrained case its rate is
  still exact for. It has to be written differently from the others: it
  terminates on time rather than on events, because `DS`'s jump chain is
  uniform over the legal moves and so carries no model information — an
  `events=` endpoint returns near-uniform draws whatever `coef` says — and it
  starts each run from a draw from the target, so what is tested is that the
  process preserves the distribution.
