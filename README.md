# `ergmgp`: Tools for Modeling ERGM Generating Processes

> **This is a fork.** It adds two features to
> [`statnet/ergmgp`](https://github.com/statnet/ergmgp), each on its own branch,
> and tracks upstream otherwise. `main` is upstream, unmodified — install from a
> branch, not from here.
>
> | branch | adds | install |
> |---|---|---|
> | [`constraints`](../../tree/constraints) | constrained sample spaces: `~edges`, `~degrees`, `~odegrees`, `~idegrees`, and dyad-level constraints | `remotes::install_github("benrosche/ergmgp@constraints")` |
> | [`rate-esp`](../../tree/rate-esp) | the above **plus** `rate.esp`/`rate.esp.cap`, a dyad-varying pacing factor | `remotes::install_github("benrosche/ergmgp@rate-esp")` |
>
> `rate-esp` is stacked on `constraints`, so it contains everything. `ergm` does
> not need reinstalling; the package compiles against the headers already in your
> library. On Windows you need Rtools, since this builds C code. Pin the ref
> explicitly in cluster scripts — these branches are rebased.
>
> **What the branches do.** A constraint restricts the sample space, which under
> a count-preserving constraint changes what an *event* is: a move has to toggle
> two or four dyads at once, which extends the process class of Butts (2023) and
> not merely its implementation. `rate.esp` multiplies a dyad's rate by
> `exp(λ·min(ESP, cap))`, which provably leaves the ERGM equilibrium untouched
> and so separates an *opportunity* effect (what a friend recommender does to a
> choice set) from a *preference* effect (a coefficient change).
>
> **Three things worth knowing before you run anything:**
>
> - **Terminate on `time=`, not `events=`.** Event times are not random times, so
>   a fixed event count biases the final state. For `DS` it is worse than biased:
>   its jump chain is uniform over the legal moves and carries no information
>   about the model at all.
> - **`DS` under `~degrees` is refused.** Its rate carries `1/|H|`, so its
>   equilibrium there is `|H|·exp(pot)` rather than the ERGM you asked for. Use
>   `LERGM` or `CI`; `DS` is exact under every other constraint.
> - **`rate.esp` needs single-dyad moves**, so it is unavailable under any
>   count-preserving constraint, and its cap defaults to 5 — uncapped, the
>   thinning bound must assume the most-connected dyad and the feature does not
>   scale past a few hundred nodes.
>
> `EGPConstraintSupport()` prints the full process-by-constraint table with the
> reason for every unsupported cell, and `NEWS.md` on either branch has the
> derivations. The branches are intended as two separate pull requests upstream.

---

`ergmgp` contains tools for modeling continuous time graph processes with equilibrium distributions in exponential family random graph model (ERGM) form (*ERGM generating processes*, or *EGPs*).  The `ergmgp` tools allow models to be specified in terms of their equilibria, using [`ergm`](https://github.com/statnet/ergm) formulas, together with additional information on the dynamic process.  A number of different processes are supported, including longitudinal ERGMs, continuum limits of temporal and separable temporal ERGMs, and competing rate stochastic actor-oriented model processes, as well as differential stability and change inhibition processes.

An overview of supported EGPs and pointers to other help pages can be obtained after loading the package with `help(ergmgp)`.  Additional information on EGPs can also be found at the reference below.

## Installing from CRAN

The easiest way to install the package is to use CRAN.  From within `R`, simply use

```
install.packages("ergmgp")
```

which will install `ergmgp` and its dependencies.  Calling `library(ergmgp)` will subsequently load the package, and away you go.

## Installing Directly from GitHub

To install from GitHub, first ensure that you have the `devtools` package installed and loaded. Then, type the following: 

```
install_github("statnet/ergmgp")
```
Alternately, cloning this repository and building/installing the package locally is another option. 

## References

Butts, Carter T.  (2023).  ``Continuous Time Graph Processes with Known ERGM Equilibria: Contextual Review, Extensions, and Synthesis.'' *Journal of Mathematical Sociology*.  DOI: 10.1080/0022250X.2023.2180001  [[arXiv version](https://arxiv.org/abs/2203.06948)]
