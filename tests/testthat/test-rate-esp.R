#Tests for rate.esp, the dyad-varying pacing factor.
#
#The claim that has to hold is that a symmetric rate modulation changes the
#dynamics but not the equilibrium.  If that fails, the feature is wrong and
#should not be used: it would silently be simulating a different ERGM than the
#one the user specified.

esp_net <- function(n = 24, nedge = 60, seed = 1){
  set.seed(seed)
  nw <- network.initialize(n, directed = FALSE)
  nw %v% "sex" <- rep(0:1, each = n/2)
  e <- combn(n, 2)
  s <- sample(ncol(e), nedge)
  add.edges(nw, e[1, s], e[2, s])
}


test_that("rate.esp is refused where the symmetry argument does not hold", {
  nw <- esp_net()
  f <- nw ~ edges + gwesp(0.3, fixed = TRUE)
  #Multi-dyad moves: ESP is not invariant, so the modulation is not symmetric.
  #This holds for the 2-toggle ~edges swap as much as for the ~degrees tetrad --
  #the dissolved edge is part of the formed dyad's neighbourhood in one state
  #and not the other -- so both are refused.
  for(con in list(~degrees, ~edges))
    expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                        constraints = con, rate.esp = 0.5, verbose = FALSE),
                 "single-dyad moves", label = deparse(con))
  #Directed: shared-partner count is ambiguous
  nd <- network.initialize(12, directed = TRUE)
  nd <- add.edges(nd, c(1,2,3,4), c(2,3,4,5))
  expect_error(simEGP(nd ~ edges, coef = -1, events = 5, process = "LERGM",
                      rate.esp = 0.5, verbose = FALSE),
               "undirected networks only")
  #Nonsense values
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      rate.esp = c(1, 2), verbose = FALSE), "single finite number")
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      rate.esp = Inf, verbose = FALSE), "single finite number")
})


test_that("rate.esp = 0 is exactly the unmodulated process", {
  nw <- esp_net()
  f <- nw ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE)
  set.seed(4); a <- simEGP(f, coef = c(-2, .6, .5), events = 200,
                           process = "LERGM", verbose = FALSE)
  set.seed(4); b <- simEGP(f, coef = c(-2, .6, .5), events = 200,
                           process = "LERGM", rate.esp = 0, verbose = FALSE)
  expect_identical(as.matrix(a), as.matrix(b))
  expect_equal(a %n% "Time", b %n% "Time")
})


test_that("rate.esp changes the dynamics", {
  #Embedded dyads are reconsidered more often, so in a fixed span of time the
  #process should take more steps.
  skip_on_cran()
  nw <- esp_net()
  f <- nw ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE)
  ev <- function(lam){
    set.seed(21)
    mean(vapply(1:25, function(i)
      simEGP(f, coef = c(-2.2, .8, .6), time = 8, process = "LERGM",
             rate.esp = lam, verbose = FALSE) %n% "Events", numeric(1)))
  }
  expect_gt(ev(0.6), ev(0) * 1.2)
})


test_that("rate.esp leaves the equilibrium where it was", {
  #The load-bearing test.  Same model, same coefficients, very different pacing:
  #the equilibrium distribution must be unchanged, because the modulation is
  #symmetric in the two states each toggle connects.
  #
  #Keep this cheap.  rate.esp rescales the clock -- the total event rate rises
  #by roughly exp(lambda * mean ESP) -- so running both arms to the same `time`
  #makes the high-lambda arm do an order of magnitude more work.  Cost is
  #superlinear in lambda and in the horizon, and an earlier version of this test
  #(lambda = 0.8, time = 60, 300 reps) ran for hours.  The claim needs both arms
  #*equilibrated*, not run long, and the high-lambda arm equilibrates faster in
  #clock time, so a short horizon is if anything more favourable to it.
  skip_on_cran()
  co <- c(-2.2, 0.8, 0.6)
  mono <- function(x) summary(x ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE))
  reps <- 150
  run <- function(lam){
    set.seed(99)
    t(vapply(seq_len(reps), function(i){
      nw <- esp_net(seed = i)
      mono(simEGP(nw ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE),
                  coef = co, time = 25, process = "LERGM",
                  rate.esp = lam, verbose = FALSE))
    }, numeric(3)))
  }
  A <- run(0); B <- run(0.5)
  se <- sqrt(apply(A, 2, var)/reps + apply(B, 2, var)/reps)
  z  <- (colMeans(A) - colMeans(B)) / se
  expect_lt(max(abs(z)), 4)
})
