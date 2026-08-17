#Tests for constrained EGP simulation.
#
#The tests are ordered from cheapest and most decisive to most expensive.  The
#two that actually earn their keep are:
#
#  - "the two engines agree on the mean inter-event time", which is what caught
#    a factor-of-two error in the tetrad move set (each move was being generated
#    twice by the sampler, doubling the total rate); and
#  - "the constrained equilibrium matches ergm::simulate", which is what would
#    catch a mistake in the multi-toggle change statistics.
#
#Both are run against ~degrees and ~edges.  Under ~edges the equilibrium check
#also monitors a degree count, which is the marginal ~degrees pins and ~edges
#frees -- so it is the one a wrong 2-toggle move set could get wrong while still
#reproducing the model terms.

undirected_net <- function(n = 14, nedge = 28, seed = 42){
  set.seed(seed)
  nw <- network.initialize(n, directed = FALSE)
  nw %v% "sex" <- rep(0:1, length.out = n)
  e <- combn(n, 2)
  s <- sample(ncol(e), nedge)
  add.edges(nw, e[1, s], e[2, s])
}

directed_net <- function(n = 16, p = 0.2, seed = 2){
  set.seed(seed)
  nw <- network.initialize(n, directed = TRUE)
  e <- which(matrix(runif(n * n), n) < p, arr.ind = TRUE)
  e <- e[e[, 1] != e[, 2], ]
  add.edges(nw, e[, 1], e[, 2])
}


test_that("unrecognised arguments are an error, not a silent no-op", {
  nw <- undirected_net()
  f <- nw ~ edges + gwesp(0.3, fixed = TRUE)
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      bogus = 1, verbose = FALSE),
               "Unused argument")
  #`constraint` (singular) partial-matches `constraints`, so it now applies the
  #constraint rather than being swallowed.
  d0 <- summary(nw ~ degree(0:10))
  s <- simEGP(f, coef = c(-2, .4), events = 50, process = "LERGM",
              constraint = ~degrees, verbose = FALSE)
  expect_identical(summary(s ~ degree(0:10)), d0)
})


test_that("unsupported process/constraint combinations are refused explicitly", {
  nw <- undirected_net()
  f <- nw ~ edges + gwesp(0.3, fixed = TRUE)
  #Continuum STERGMs cannot take multi-toggle moves at all.  (CDCSTERGM takes a
  #single formula with a list-valued coef, so it reaches the constraint check;
  #CSTERGM proper needs a pair of formulas and would fail earlier, on the model.)
  expect_error(simEGP(f, coef = list(formation = c(-2, .4), dissolution = -1),
                      events = 5, process = "CDCSTERGM",
                      constraints = ~degrees, verbose = FALSE),
               "separability")
  #A constraint that freezes every dyad is a specification error, not a
  #simulation that silently never moves
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      constraints = ~blocks(attr = "sex", levels2 = TRUE),
                      verbose = FALSE),
               "no free dyads")
  #Constraints we cannot turn into a move set.  (~degreedist fixes the degree
  #*distribution* rather than the degree sequence, so no fixed-shape move set
  #preserves it; it exposes no free-dyad set either.)
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      constraints = ~degreedist, verbose = FALSE),
               "cannot simulate under the constraint")
  #Only one count-preserving constraint at a time
  expect_error(simEGP(directed_net() ~ edges + mutual, coef = c(-1, .5), events = 5,
                      process = "LERGM", constraints = ~odegrees + idegrees,
                      verbose = FALSE),
               "Only one count-preserving constraint")
  #Thinning needs a bounded rate
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "CTERGM",
                      constraints = ~degrees, engine = "thinning", verbose = FALSE),
               "not available for process")
  #DS needs the move count, so it cannot be thinned under a count-preserving
  #constraint
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "DS",
                      constraints = ~edges, engine = "thinning", verbose = FALSE),
               "not available for DS")
  #...and under ~degrees the move count is not merely unknown but state
  #dependent, which moves DS's equilibrium off the requested ERGM, so the
  #combination is refused whatever engine is asked for
  for(eng in c("auto", "thinning", "enumeration"))
    expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "DS",
                        constraints = ~degrees, engine = eng, verbose = FALSE),
                 "DS cannot be simulated under ~degrees", fixed = TRUE)
  #The refusal is specific to the process and to the constraint, not general to
  #either: DS runs under the other count-preserving constraints, and the other
  #processes run under ~degrees.
  expect_s3_class(simEGP(f, coef = c(-2, .4), events = 5, process = "DS",
                         constraints = ~edges, verbose = FALSE), "network")
  expect_s3_class(simEGP(f, coef = c(-2, .4), events = 5, process = "CI",
                         constraints = ~degrees, verbose = FALSE), "network")
  #Directedness mismatches
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      constraints = ~odegrees, verbose = FALSE))
})


test_that("degree-preserving constraints preserve exactly what they claim", {
  nw <- undirected_net()
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  d0 <- summary(nw ~ degree(0:12))
  for(eng in c("thinning", "enumeration")){
    s <- simEGP(f, coef = c(0.6, 0.5), events = 200, process = "LERGM",
                constraints = ~degrees, engine = eng, verbose = FALSE)
    expect_identical(summary(s ~ degree(0:12)), d0, info = eng)
    expect_identical(summary(s ~ edges), summary(nw ~ edges), info = eng)
    expect_equal(s %n% "Engine", eng)
  }
  nd <- directed_net()
  fd <- nd ~ edges + mutual
  od0 <- summary(nd ~ odegree(0:10)); id0 <- summary(nd ~ idegree(0:10))
  s <- simEGP(fd, coef = c(-1.5, .5), events = 300, process = "LERGM",
              constraints = ~odegrees, verbose = FALSE)
  expect_identical(summary(s ~ odegree(0:10)), od0)   #Held fixed
  expect_false(identical(summary(s ~ idegree(0:10)), id0))  #Free to move
  s <- simEGP(fd, coef = c(-1.5, .5), events = 300, process = "LERGM",
              constraints = ~idegrees, verbose = FALSE)
  expect_identical(summary(s ~ idegree(0:10)), id0)
})


test_that("~edges fixes the edge count and frees the degree sequence", {
  #This is the whole point of the constraint: it sits between ~degrees (which
  #forbids any change to the degree distribution) and the unconstrained process
  #(which lets volume and composition move together).  So the test has to check
  #both halves -- that the count is held, *and* that the degree sequence is not.
  nw <- undirected_net()
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  d0 <- summary(nw ~ degree(0:12)); e0 <- summary(nw ~ edges)
  for(eng in c("thinning", "enumeration")){
    s <- simEGP(f, coef = c(0.6, 0.5), events = 200, process = "LERGM",
                constraints = ~edges, engine = eng, verbose = FALSE)
    expect_identical(summary(s ~ edges), e0, info = eng)
    expect_false(identical(summary(s ~ degree(0:12)), d0), info = eng)
    expect_equal(s %n% "Engine", eng)
    expect_equal(s %n% "Constraint", "~edges")
  }
  #Unlike the degree constraints, ~edges applies to directed networks too
  nd <- directed_net()
  s <- simEGP(nd ~ edges + mutual, coef = c(-1.5, .5), events = 300,
              process = "LERGM", constraints = ~edges, verbose = FALSE)
  expect_identical(summary(s ~ edges), summary(nd ~ edges))
  expect_false(identical(summary(s ~ odegree(0:10)), summary(nd ~ odegree(0:10))))
  #ergm resolves the implication for us: ~degrees fixes every degree and so
  #fixes their sum, and ergm_conlist() drops the redundant ~edges.  Adding it
  #should therefore be a no-op rather than a conflict.
  s <- simEGP(f, coef = c(0.6, 0.5), events = 50, process = "LERGM",
              constraints = ~edges + degrees, verbose = FALSE)
  expect_equal(s %n% "Constraint", "~degrees")
  expect_identical(summary(s ~ degree(0:12)), d0)
})


test_that("under ~edges every event is one removal and one addition", {
  #The move is a 2-toggle swap, so the event history must come in balanced
  #pairs sharing a timestamp.  extract_stats()-style analyses depend on this.
  nw <- undirected_net()
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  s <- simEGP(f, coef = c(0.6, 0.5), events = 40, process = "LERGM",
              constraints = ~edges, return.history = TRUE, verbose = FALSE)
  eh <- s %n% "EventHistory"
  expect_equal(nrow(eh), 2 * (s %n% "Events"))
  expect_equal(sum(eh[, "Onset"] == 1), sum(eh[, "Onset"] == 0))
  expect_true(all(tapply(eh[, "Onset"], eh[, "Time"],
                         function(x) length(x) == 2L && sum(x) == 1)))
})


test_that("dyad-level constraints toggle only free dyads", {
  set.seed(9)
  nw <- undirected_net()
  nw %v% "grp" <- rep(1:2, length.out = network.size(nw))
  f <- nw ~ edges + gwesp(0.3, fixed = TRUE)
  #levels2 names the block pairs to *forbid*; 2 is the cross-group pair, so this
  #leaves only within-group dyads free (42 of 91 here).
  s <- simEGP(f, coef = c(-1.5, .4), events = 300, process = "LERGM",
              constraints = ~blocks(attr = "grp", levels2 = 2), verbose = FALSE)
  #No cross-group edge may have been created or destroyed
  a <- as.matrix(nw); b <- as.matrix(s)
  g <- nw %v% "grp"
  cross <- outer(g, g, "!=")
  expect_identical(a[cross], b[cross])
})


test_that("the two engines agree on the mean inter-event time", {
  #This is the sharpest available check on the rate algebra: the embedded jump
  #chain is insensitive to a constant factor in the rates, but the waiting time
  #is not.  A move set that generates some moves twice shows up here and
  #nowhere else.
  skip_on_cran()
  nw <- undirected_net()
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  reps <- 120; ev <- 100
  mit <- function(eng, con){
    set.seed(3)
    mean(vapply(1:reps, function(i)
      simEGP(f, coef = c(0.6, 0.5), events = ev, process = "LERGM",
             constraints = con, engine = eng, verbose = FALSE) %n% "Time",
      numeric(1))) / ev
  }
  for(con in list(~degrees, ~edges)){
    a <- mit("thinning", con); b <- mit("enumeration", con)
    expect_equal(a / b, 1, tolerance = 0.06, label = deparse(con))
  }
})


test_that("the constrained equilibrium matches ergm::simulate", {
  #Both target the same constrained ERGM, so their equilibrium means must agree
  #up to Monte Carlo error.  This is what catches errors in the multi-toggle
  #change statistics, which is where gwesp's auxiliary storage could bite.
  skip_on_cran()
  nw <- undirected_net()
  co <- c(0.6, 0.5)
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  reps <- 250
  #For ~edges we also monitor a degree count.  That is precisely the marginal
  #~degrees pins and ~edges frees, so it is the one a wrong 2-toggle move set
  #would get wrong while still reproducing the model terms.  It is degenerate
  #under ~degrees, so it cannot be monitored there.
  cases <- list(list(con = ~degrees, mon = "gwesp(0.3, fixed = TRUE) + nodematch(\"sex\")"),
                list(con = ~edges,   mon = "gwesp(0.3, fixed = TRUE) + nodematch(\"sex\") + degree(2)"))
  for(cs in cases){
    mono <- function(x) summary(as.formula(paste("x ~", cs$mon)))
    k <- length(mono(nw))
    set.seed(7)
    A <- t(vapply(1:reps, function(i)
      mono(simEGP(f, coef = co, events = 4000, process = "LERGM",
                  constraints = cs$con, verbose = FALSE)), numeric(k)))
    set.seed(11)
    B <- t(vapply(1:reps, function(i)
      mono(suppressMessages(
        simulate(f, coef = co, constraints = cs$con, nsim = 1,
                 control = control.simulate.formula(MCMC.burnin = 20000),
                 output = "network"))), numeric(k)))
    se <- sqrt(apply(A, 2, var)/reps + apply(B, 2, var)/reps)
    z <- (colMeans(A) - colMeans(B)) / se
    expect_lt(max(abs(z)), 4, label = deparse(cs$con))
  }
})


test_that("DS keeps the constrained ERGM where the move count is fixed", {
  #The other half of refusing DS under ~degrees: the claim that it is exact
  #wherever the constraint pins |H|, here ~edges, where |H| = E(D-E).
  #
  #Two things make this test unlike the one above.  It must terminate on TIME:
  #DS's rate does not depend on which move is taken, so its jump chain is
  #uniform over the legal moves and carries no model information at all, and an
  #events= endpoint would return draws that are near-uniform over the reachable
  #states however the coefficients are set.  The model enters only through the
  #holding times.  And it starts each run *from* the target rather than from
  #nw, so what is being tested is whether the process preserves the
  #distribution -- which is the property at issue -- with no burn-in to argue
  #about.  DS's exit rate is A exp(-pot), so the pacing constant has to undo
  #exp(pot) or nothing happens in any reachable amount of time.
  skip_on_cran()
  nw <- undirected_net()
  co <- c(0.6, 0.5)
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  mono <- function(x)
    summary(x ~ gwesp(0.3, fixed = TRUE) + nodematch("sex") + degree(2))
  target <- function(){
    suppressMessages(simulate(f, coef = co, constraints = ~edges, nsim = 1,
                              control = control.simulate.formula(MCMC.burnin = 20000),
                              output = "network"))
  }
  reps <- 200
  set.seed(7)
  A <- t(vapply(1:reps, function(i){
    d <- target()
    rf <- exp(sum(co * summary(d ~ gwesp(0.3, fixed = TRUE) + nodematch("sex"))))
    s <- simEGP(d ~ gwesp(0.3, fixed = TRUE) + nodematch("sex"), coef = co,
                time = 100, rate.factor = rf, process = "DS",
                constraints = ~edges, verbose = FALSE)
    c(mono(s), s %n% "Events")
  }, numeric(4)))
  set.seed(11)
  B <- t(vapply(1:reps, function(i) mono(target()), numeric(3)))
  se <- sqrt(apply(A[, 1:3], 2, var)/reps + apply(B, 2, var)/reps)
  expect_lt(max(abs((colMeans(A[, 1:3]) - colMeans(B)) / se)), 4)
  #Guard against the run being vacuous, which is how this test would silently
  #stop testing anything -- but in aggregate, because no per-run floor is
  #achievable here.  The exit rate is A exp(-pot) and the potential spans about
  #ten nats across the target, so the event count per unit time spans four
  #orders of magnitude; pacing hard enough to put a floor under the slowest
  #draw makes the fastest one run for millions of events.  A run that barely
  #moves costs power, not validity: its endpoint is still a draw from the
  #target, which is exactly the null being tested.
  expect_gt(median(A[, 4]), 100)
})
