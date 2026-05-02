# Univariate prior on model parameter

A `Prior` object represents a prior distribution on the single model
parameter of a
[`DataDistribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
class object. Together a prior and data-distribution specify the class
of the joint distribution of the test statisic, X, and its parameter,
theta. Currently, adoptr only allows simple models with a single
parameter. Implementations for
[PointMassPrior](https://optad.github.io/adoptr/reference/PointMassPrior-class.md)
and
[ContinuousPrior](https://optad.github.io/adoptr/reference/ContinuousPrior-class.md)
are available.

## Details

For an example on working with priors, see
[here](https://optad.github.io/adoptr/articles/working-with-priors.html).

## See also

For the available methods, see
[`bounds`](https://optad.github.io/adoptr/reference/bounds.md),
[`expectation`](https://optad.github.io/adoptr/reference/expectation.md),
[`condition`](https://optad.github.io/adoptr/reference/condition.md),
[`predictive_pdf`](https://optad.github.io/adoptr/reference/predictive_pdf.md),
[`predictive_cdf`](https://optad.github.io/adoptr/reference/predictive_cdf.md),
[`posterior`](https://optad.github.io/adoptr/reference/posterior.md)

## Examples

``` r
disc_prior <- PointMassPrior(c(0.1, 0.25), c(0.4, 0.6))

cont_prior <- ContinuousPrior(
  pdf     = function(x) dnorm(x, mean = 0.3, sd = 0.2),
  support = c(-2, 3)
)

```
