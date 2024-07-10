context("Chi-Squared distribution")

test_that('constructors work', {
  
  expect_true(
    ChiSquared(df = 10)@multiplier == 1
  )
  
  expect_error(
    ChiSquared(df = -1)
  )
  
  expect_error(
    ChiSquared(df = 2.3)
  )
  
  expect_true(
    Pearson2xK(5)@df == 4
  )
  
  expect_true(
    Pearson2xK(5)@multiplier == 1 / (Pearson2xK(5)@df + 1)
  )
  
  expect_error(
    Pearson2xK(-3)
  )
  
  expect_error(
    Pearson2xK(1.8)
  )
  
  expect_equal(
    ZSquared()@multiplier, ZSquared()@df + 1
  )
  
  expect_equal(
    ZSquared(FALSE)@multiplier, ZSquared()@df
  )
})

test_that('thetas and ncps are correctly computed', {
  
  # theoretical ncp and calculated ncp are equal
  # Pearson2xK
  p_vector <- c(0.8, 0.15, 0.3, 0.96)
  n <- 52
  dist <- Pearson2xK(length(p_vector))
  real_ncp <- n / length(p_vector) * (sum((p_vector - mean(p_vector))^2)) /
    (mean(p_vector) * (1 - mean(p_vector)))
  ncp_calc <- n / dist@multiplier * get_tau_Pearson2xK(p_vector) * 1 / length(p_vector)
  
  expect_equal(ncp_calc, real_ncp)
  
  #ZSquared
  dist2 <- ZSquared(TRUE)
  mu <- 0.4
  sigma <- 1.2
  expect_equal(
    get_tau_ZSquared(mu, sigma) * n / dist2@multiplier,
    (sqrt(n / 2) * mu/sigma)^2
  )
  
  dist3 <- ZSquared(FALSE)
  expect_equal(
    get_tau_ZSquared(mu, sigma) * n / dist3@multiplier,
    (sqrt(n) * mu/sigma)^2
  )
})

test_that('pdf is defined correctly', {
  
  dist <- Pearson2xK(3)
  x       <- seq(-3, 3, by = .1)
  n       <- 22
  p_vec   <- c(0.5, 0.6, 0.7)
  theta   <- get_tau_Pearson2xK(p_vec)
  expect_equal(
    probability_density_function(dist, x, n, theta),
    stats::dchisq(x, df = 2, ncp = n * (sum((p_vec - mean(p_vec))^2)) /
                    (mean(p_vec) * (1 - mean(p_vec)))),
    tolerance = 1e-6, scale = 1)
})

test_that('cdf is defined correctly', {
  dist <- Pearson2xK(3)
  x       <- seq(-3, 3, by = .1)
  n       <- 22
  p_vec   <- c(0.5, 0.6, 0.7)
  theta   <- get_tau_Pearson2xK(p_vec)
  expect_equal(
    cumulative_distribution_function(dist, x, n, theta),
    stats::pchisq(x, df = 2, ncp = n * (sum((p_vec - mean(p_vec))^2)) /
                    (mean(p_vec) * (1 - mean(p_vec)))),
    tolerance = 1e-6, scale = 1)
})

test_that('quantile is defined correctly', {
  
  dist <- Pearson2xK(3)
  x       <- seq(0.1, 1, by = .01)
  n       <- 22
  p_vec   <- c(0.5, 0.6, 0.7)
  theta   <- get_tau_Pearson2xK(p_vec)
  expect_equal(
    quantile(dist, x, n, theta),
    stats::qchisq(x, df = 2, ncp = n * (sum((p_vec - mean(p_vec))^2)) /
                    (mean(p_vec) * (1 - mean(p_vec)))),
    tolerance = 1e-6, scale = 1)
})

test_that('simulate respects seed', {
  
  expect_equal(
    simulate(Pearson2xK(3), 10, 10, 0.2, seed = 42),
    simulate(Pearson2xK(3), 10, 10, 0.2, seed = 42),
    tolerance = 1e-6, scale = 1)
  
  set.seed(42)
  
  expect_true(
    all(simulate(Pearson2xK(4), 10, 12, 1.1) != simulate(Pearson2xK(4), 10, 12, 1.1)))
})

test_that("show method", {
  
  expect_equal(
    capture.output(show(Pearson2xK(4))),
    "Pearson2xK<df=3> "
  )
  
  expect_equal(
    capture.output(show(ZSquared())),
    "ZSquared<df=1> "
  )
  
  expect_equal(
    capture.output(show(ChiSquared(5))),
    "ChiSquared<df=5> "
  )
})