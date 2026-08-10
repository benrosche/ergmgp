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
  #Constraints we cannot turn into a move set
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      constraints = ~edges, verbose = FALSE),
               "cannot simulate under the constraint")
  #Thinning needs a bounded rate
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "CTERGM",
                      constraints = ~degrees, engine = "thinning", verbose = FALSE),
               "not available for process")
  #DS needs the move count, so it cannot be thinned under a degree constraint
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "DS",
                      constraints = ~degrees, engine = "thinning", verbose = FALSE),
               "not available for DS")
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
  mit <- function(eng){
    set.seed(3)
    mean(vapply(1:reps, function(i)
      simEGP(f, coef = c(0.6, 0.5), events = ev, process = "LERGM",
             constraints = ~degrees, engine = eng, verbose = FALSE) %n% "Time",
      numeric(1))) / ev
  }
  a <- mit("thinning"); b <- mit("enumeration")
  expect_equal(a / b, 1, tolerance = 0.06)
})


test_that("the constrained equilibrium matches ergm::simulate", {
  #Both target the same constrained ERGM, so their equilibrium means must agree
  #up to Monte Carlo error.  This is what catches errors in the multi-toggle
  #change statistics, which is where gwesp's auxiliary storage could bite.
  skip_on_cran()
  nw <- undirected_net()
  co <- c(0.6, 0.5)
  f <- nw ~ gwesp(0.3, fixed = TRUE) + nodematch("sex")
  mono <- function(x) summary(x ~ gwesp(0.3, fixed = TRUE) + nodematch("sex"))
  reps <- 250
  set.seed(7)
  A <- t(vapply(1:reps, function(i)
    mono(simEGP(f, coef = co, events = 4000, process = "LERGM",
                constraints = ~degrees, verbose = FALSE)), numeric(2)))
  set.seed(11)
  B <- t(vapply(1:reps, function(i)
    mono(simulate(f, coef = co, constraints = ~degrees, nsim = 1,
                  control = control.simulate.formula(MCMC.burnin = 20000),
                  output = "network")), numeric(2)))
  se <- sqrt(apply(A, 2, var)/reps + apply(B, 2, var)/reps)
  z <- (colMeans(A) - colMeans(B)) / se
  expect_lt(max(abs(z)), 4)
})
