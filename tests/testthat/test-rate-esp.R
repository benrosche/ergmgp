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

#A network with one high-degree node, for the tests about the rate bound: the
#uncapped bound is driven by the maximum degree, so a hub is what separates the
#capped and uncapped cases.  Keep hubdeg modest -- the uncapped arm's acceptance
#is about exp(-rate.esp * hubdeg), and it has to stay slow-but-finite.
hub_net <- function(n = 60, hubdeg = 20, nedge = 90, seed = 7){
  set.seed(seed)
  nw <- network.initialize(n, directed = FALSE)
  nw <- add.edges(nw, rep(1, hubdeg), 2:(hubdeg + 1))
  e <- combn(2:n, 2)
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
  #rate.esp.cap=Inf is given explicitly so that this test keeps meaning what it
  #was written to mean once the cap gained a finite default; the capped variant
  #is covered separately below.
  ev <- function(lam){
    set.seed(21)
    mean(vapply(1:25, function(i)
      simEGP(f, coef = c(-2.2, .8, .6), time = 8, process = "LERGM",
             rate.esp = lam, rate.esp.cap = Inf, verbose = FALSE) %n% "Events", numeric(1)))
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
                  rate.esp = lam, rate.esp.cap = Inf, verbose = FALSE))
    }, numeric(3)))
  }
  A <- run(0); B <- run(0.5)
  se <- sqrt(apply(A, 2, var)/reps + apply(B, 2, var)/reps)
  z  <- (colMeans(A) - colMeans(B)) / se
  expect_lt(max(abs(z)), 4)
})


test_that("rate.esp.cap is validated, including when the pacing is off", {
  nw <- esp_net()
  f <- nw ~ edges + gwesp(0.3, fixed = TRUE)
  bad <- function(cap, lam = 0.5)
    expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                        rate.esp = lam, rate.esp.cap = cap, verbose = FALSE),
                 "non-negative")
  #A negative cap is not merely meaningless: min(ESP,cap) would be constant and
  #negative, putting the true rate above the bound the thinning engine assumes.
  bad(-1); bad(c(1, 2)); bad(NA_real_); bad(NA)
  #Checked before rate.esp==0 short-circuits, so a bad cap is never accepted
  #just because the modulation happens to be switched off.
  bad(-1, lam = 0)
  #rate.esp itself is still finite-only, and the two messages stay distinct
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      rate.esp = Inf, rate.esp.cap = 5, verbose = FALSE),
               "single finite number")
  #The correctness restrictions still fire when a cap is in force
  expect_error(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                      constraints = ~edges, rate.esp = .5, rate.esp.cap = 3,
                      verbose = FALSE), "single-dyad moves")
})


test_that("rate.esp.cap = 0 is exactly the unmodulated process", {
  #min(ESP,0) = 0, so the modulation is identically exp(0) = 1.  This has to hold
  #bit-for-bit, not just in distribution: it pins that the cap enters the rate
  #and the bound as an exact 0.0 rather than as something merely small.  Compare
  #at a fixed engine, since the engine rule legitimately differs between arms.
  nw <- esp_net()
  f <- nw ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE)
  for(eng in c("enumeration", "thinning")){
    set.seed(4); a <- simEGP(f, coef = c(-2, .6, .5), events = 200, process = "LERGM",
                             engine = eng, verbose = FALSE)
    set.seed(4); b <- simEGP(f, coef = c(-2, .6, .5), events = 200, process = "LERGM",
                             engine = eng, rate.esp = 0.7, rate.esp.cap = 0,
                             verbose = FALSE)
    expect_identical(as.matrix(a), as.matrix(b), info = eng)
    expect_equal(a %n% "Time", b %n% "Time")
    expect_null(b %n% "RateESP")   #Normalised away, so the record is honest
  }
})


test_that("a cap above the achievable ESP changes nothing", {
  #esp_net() has 24 nodes, so ESP <= 22 and a cap of 100 can never bind.  Fixing
  #the engine is necessary: under thinning the two arms have different bounds
  #(resp*100 versus resp*max degree), so their RNG streams diverge even though
  #both are exact.
  nw <- esp_net()
  f <- nw ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE)
  set.seed(8); a <- simEGP(f, coef = c(-2, .6, .5), events = 200, process = "LERGM",
                           engine = "enumeration", rate.esp = 0.5,
                           rate.esp.cap = Inf, verbose = FALSE)
  set.seed(8); b <- simEGP(f, coef = c(-2, .6, .5), events = 200, process = "LERGM",
                           engine = "enumeration", rate.esp = 0.5,
                           rate.esp.cap = 100, verbose = FALSE)
  expect_identical(as.matrix(a), as.matrix(b))
  expect_equal(a %n% "Time", b %n% "Time")
})


test_that("the cap decides the engine, however `auto` is spelt", {
  #The rule: with a pacing factor in force, "auto" thins iff the cap is finite,
  #because only then is the bound tight.  This also guards a bug that used to
  #live in simEGP(): the fallback to enumeration tested `identical(engine,
  #c("auto","thinning","enumeration"))`, i.e. the unevaluated default, so it did
  #not fire for an explicit engine="auto" nor for anything arriving from
  #simEGPTraj() -- which under a dyad-level constraint silently used the loose
  #bound.  The rule now lives in EGP_engine() and applies to all three spellings.
  nw <- esp_net()
  f <- nw ~ edges + gwesp(0.3, fixed = TRUE)
  eng <- function(...) suppressWarnings(
    simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM", verbose = FALSE,
           ...) %n% "Engine")
  expect_equal(eng(rate.esp = 0), "enumeration")
  expect_equal(eng(rate.esp = 0.5, rate.esp.cap = Inf), "enumeration")
  expect_equal(eng(rate.esp = 0.5, rate.esp.cap = Inf, engine = "auto"), "enumeration")
  expect_equal(eng(rate.esp = 0.5, rate.esp.cap = 4), "thinning")
  expect_equal(eng(rate.esp = 0.5), "thinning")            #The default cap is finite
  #The uncapped thinning path is reachable on request, and says why it is a bad idea
  expect_warning(simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
                        engine = "thinning", rate.esp = 0.5, rate.esp.cap = Inf,
                        verbose = FALSE),
                 "most-connected dyad")
  #The run records the model it was actually simulated under, so a verbose=FALSE
  #caller is not left guessing which cap applied
  s <- simEGP(f, coef = c(-2, .4), events = 5, process = "LERGM",
              rate.esp = 0.5, verbose = FALSE)
  expect_equal(s %n% "RateESP", 0.5)
  expect_equal(s %n% "RateESPCap", 5)
  #simEGPTraj must forward the cap and validate up front rather than inside a
  #worker; forgetting the former would silently run uncapped trajectories.
  tr <- simEGPTraj(f, coef = c(-2, .4), events = 5, process = "LERGM",
                   rate.esp = 0.5, rate.esp.cap = 4, verbose = FALSE)
  expect_equal(tr[[2]] %n% "Engine", "thinning")
  expect_equal(tr[[2]] %n% "RateESPCap", 4)
  expect_error(simEGPTraj(f, coef = c(-2, .4), events = 5, process = "LERGM",
                          constraints = ~edges, rate.esp = 0.5, verbose = FALSE),
               "single-dyad moves")
})


test_that("a finite cap makes the thinning engine viable", {
  #The point of the cap.  Uncapped, the bound has to assume the most-connected
  #dyad, so acceptance is about exp(-rate.esp * max degree) -- here exp(-6), and
  #on a real network exp(-100), which is zero.  Capped, the bound is exactly
  #A*exp(rate.esp*cap) and acceptance is a workable fraction.
  skip_on_cran()
  f <- hub_net() ~ edges + gwesp(0.25, fixed = TRUE)
  acc <- function(cap){
    set.seed(11)
    suppressWarnings(
      simEGP(f, coef = c(-3, 0.3), events = 20, process = "LERGM",
             engine = "thinning", rate.esp = 0.3, rate.esp.cap = cap,
             verbose = FALSE)) %n% "Acceptance"
  }
  a_cap <- acc(3); a_unc <- acc(Inf)
  expect_gt(a_cap, 0.01)
  expect_gt(a_cap, 20 * a_unc)
})


test_that("a capped rate.esp leaves the equilibrium where it was", {
  #Capping keeps min(ESP,cap) a function of the dyad's neighbourhood alone, so
  #it is still unchanged by toggling the dyad and the detailed-balance argument
  #is untouched.  This is the same test as the uncapped one above, and it is
  #load-bearing for the same reason: if it fails, a capped run is silently
  #simulating a different ERGM than the one requested.
  #
  #cap = 2 has to bind for this to be a test of anything: esp_net() has density
  #60/276, so a typical dyad's ESP is about 1 with a tail to 5 or 6.
  skip_on_cran()
  co <- c(-2.2, 0.8, 0.6)
  mono <- function(x) summary(x ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE))
  reps <- 150
  run <- function(lam, cap){
    set.seed(99)
    t(vapply(seq_len(reps), function(i){
      nw <- esp_net(seed = i)
      mono(simEGP(nw ~ edges + nodematch("sex") + gwesp(0.3, fixed = TRUE),
                  coef = co, time = 25, process = "LERGM", rate.esp = lam,
                  rate.esp.cap = cap, verbose = FALSE))
    }, numeric(3)))
  }
  A <- run(0, Inf); B <- run(0.5, 2)
  se <- sqrt(apply(A, 2, var)/reps + apply(B, 2, var)/reps)
  expect_lt(max(abs((colMeans(A) - colMeans(B)) / se)), 4)
})
