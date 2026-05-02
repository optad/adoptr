# Predictive CDF

`predictive_cdf()` evaluates the predictive CDF of the model specified
by a
[`DataDistribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
`dist` and
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md) at
the given stage-one outcome.

## Usage

``` r
predictive_cdf(dist, prior, x1, n1, ...)

# S4 method for class 'DataDistribution,ContinuousPrior,numeric'
predictive_cdf(
  dist,
  prior,
  x1,
  n1,
  k = 10 * (prior@support[2] - prior@support[1]) + 1,
  ...
)

# S4 method for class 'DataDistribution,PointMassPrior,numeric'
predictive_cdf(dist, prior, x1, n1, ...)
```

## Arguments

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- prior:

  a [`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)
  object

- x1:

  stage-one test statistic

- n1:

  stage-one sample size

- ...:

  further optional arguments

- k:

  number of pivots for crude integral approximation

## Value

`numeric`, value of the predictive CDF

## Examples

``` r
tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
predictive_cdf(Normal(), tmp, 2, 20)
#> [1] 0.8455245

predictive_cdf(Normal(), PointMassPrior(.0, 1), 0, 20) # .5
#> [1] 0.5
```
