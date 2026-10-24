# adoptr (development version)

## Bug fixes

* `Binomial`: the null standard deviation of the test statistic is now computed per `theta`. Before, it was taken from the first element only, which gave wrong results for priors with several support points. Results for such priors change.
* `NestedModels`/`ANOVA`: pdf and cdf are now vectorized in `x`, `n` and `theta`. Before, a scalar `x` with several `theta` values (e.g. in `posterior()`) gave wrong results.
* `Student`, `NestedModels`/`ANOVA`: degrees of freedom are computed element-wise, so `simulate()` works with a vector of sample sizes.
* `summary()` for `TwoStageDesign` no longer passes the scores on to `evaluate()`, which made `summary(design, "Power" = pow)` fail.
* `n2()` for `OneStageDesign` now returns one value per `x1`, and `n2(<design>, c(-Inf, Inf))` returns `0` instead of `NaN`.
* `SurvivalDesign()` and design conversions keep the tunability status (`make_fixed()`).

## Performance

* `ConditionalPower` and `Power` with a `PointMassPrior` are computed vectorized over all stage-one points (one pdf/cdf call instead of one `posterior()`/`expectation()` per pivot). This speeds up `minimize()`.
* `n2()` and `c2()` build the spline once per call. `MaximumSampleSize` and the `plot()`/`summary()` methods use vectorized calls.
* `expand.grid()` is replaced by a lightweight grid in the predictive pdf/cdf. `minimize()` builds the design object once per parameter vector.

## Other changes

* Stricter input checks: length of `params` in `update()`, unknown names in `make_fixed()`/`make_tunable()`, scalar `n2_pivots` in `GroupSequentialDesign`, `info_ratio` in (0, 1) in `get_initial_design()`, `p_outer > p_inner` in `NestedModels`, at least two groups in `ANOVA` and `Pearson2xK`.
* Clearer errors for `posterior()`/`condition()` on a `PointMassPrior` when the observation has zero likelihood or the interval is empty. `ContinuousPrior(tighten_support = TRUE)` always terminates.
* `minimize()` also warns on negative `nloptr` status codes.
* Custom `DataDistribution` classes must be vectorized in `x`, `n` and `theta`.
* Documentation updates (distributions, scores, rounding, `bounds()` example).

# adoptr 1.1.2

* fix documentation problem according to CRAN checks
* update maintainer e-mail

# adoptr 1.1.1

* fix issue from new glue update

# adoptr 1.1.0

* include survival endpoints
* add F-distribution and chi-squared distribution
* check if constraints are fulfilled
* update `get_initial_design`
* vignettes for binomial and time-to-event endpoints

# adoptr 1.0.1

* updated maintainer
* updated roxygen version
* fulfill current CRAN requirements

# adoptr 1.0.0

* added references to article in journal of statistical software

# adoptr 0.4.3

* corrected class diagram in vignette

# adoptr 0.4.2

* minor bug/print fixes
* Student-t data distribution
* switched to GitHub Actions for CI/CD

# adoptr 0.4.0

* included JSS article as vignette

# adoptr 0.3.2

* bugfix for binomial pdf/cdf

# adoptr 0.3.1

* added support for binomial endpoint
* initial design creation more convenient
* extended summary() function for two stage designs

# adoptr 0.3.0

* bugfix for pretty printing of optimization results
* streamlined code for internal evaluation of integral scores during optimization; 
substantial speedup!

# adoptr 0.2.3

* reworked the class printing system to be more informative
* minimize() now also returns an S3 class object to allow pretty-printing



# adoptr 0.2.2

* fixed references to SIM paper



# adoptr 0.2.1

* added references to SIM paper in all relevant places



# adoptr 0.2.0

* new feature: composite scores allows more generic expressions than just 
    affine score combinations, cf. composite()
* affine scores (s1 + 2*s2) are no longer supported, use composite instead
* more consistent class system: conditional scores no longer need a specification
    of distributions by default (no need for conditional sampel size e.g.).
    Instead, expected() now requires explicit specification of the data and
    prior distribution to integrate with.
* Vignettes updated
* fixed broked tests due to updated rpact



# adoptr 0.1.1

* extended Description field in DESCRIPTION to full paragraph
* provided examples for all user facing functions
* revision of docs



# adoptr 0.1.0

* initial release
