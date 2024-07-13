context("F-distribution")

test_that("Constructors work", {
  
  expect_true(
    NestedModels(4, 5)@p_inner == 4 &
      NestedModels(4, 5)@p_outer ==5
  )
  
  expect_true(
    ANOVA(6)@p_inner == 1 &
      ANOVA(6)@p_outer == 6
  )
  
  expect_error(
    NestedModels(2.3, 1)
  )
  
  expect_error(
    NestedModels(2, -1)
  )
  
  expect_error(
    ANOVA(3.4)
  )
  
  expect_error(
    ANOVA(-6)
  )
})

test_that("thetas and ncps are correctly computed", {
  
  
  # ANOVA
  p_vec <- c(0.8, 0.15, 0.3, 0.96)
  n <- 52
  sigma <- 2
  dist <- ANOVA(length(p_vec))
  
  real_ncp <- n * (sum((p_vec - mean(p_vec))^2)) / sigma^2
  ncp_calc <- (n * (dist@p_outer - dist@p_inner + 1)) * get_tau_ANOVA(p_vec, sigma)
  
  expect_equal(
    real_ncp, ncp_calc
  )
  
})

test_that("pdf is defined correctly", {
  
  dist <- ANOVA(3)
  x <- seq(0.1, 5, by = 0.1)
  n <- 42
  p_vec <- c(0.3, 0.4, 0.5)
  theta <- get_tau_ANOVA(p_vec)
  real_ncp <- n * (sum((p_vec - mean(p_vec))^2))
  expect_equal(
    probability_density_function(dist, x, n, theta),
    stats::df(x, 2, 3 * n - 3, ncp = real_ncp)
  )
  
})

test_that("cdf is defined correctly", {
  
  dist <- ANOVA(3)
  x <- seq(0.1, 5, by = 0.1)
  n <- 32
  p_vec <- c(0.2, 0.4, 0.5)
  theta <- get_tau_ANOVA(p_vec)
  real_ncp <- n * (sum((p_vec - mean(p_vec))^2))
  expect_equal(
    cumulative_distribution_function(dist, x, n, theta),
    stats::pf(x, 2, 3 * n - 3, ncp = real_ncp)
  )
})

test_that("quantile is defined correctly", {
  
  dist <- ANOVA(4)
  x <- seq(0.1, 1, by = 0.01)
  n <- 22
  p_vec <- c(0.2, 0.4, 0.7, 0.36)
  theta <- get_tau_ANOVA(p_vec)
  real_ncp <- n * (sum((p_vec - mean(p_vec))^2))
  expect_equal(
    quantile(dist, x, n, theta),
    stats::qf(x, 3, 4 * n - 4, ncp = real_ncp)
  )
})

test_that("simulate respects seed", {
  
  expect_equal(
    simulate(ANOVA(3), 10, 10, 0.2, seed = 42),
    simulate(ANOVA(3), 10, 10, 0.2, seed = 42),
    tolerance = 1e-6, scale = 1)
  
  set.seed(42)
  
  expect_true(
    all(simulate(ANOVA(4), 10, 12, 1.1) != simulate(ANOVA(4), 10, 12, 1.1)))
})

test_that("show_method", {
  
  expect_equal(
    capture.output(show(NestedModels(4, 5))),
    "NestedModels<p_inner=4,p_outer=5> "
  )
  
  expect_equal(
    capture.output(show(ANOVA(4))),
    "ANOVA<n_groups=4> "
  )
})