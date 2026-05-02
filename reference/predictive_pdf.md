# Predictive PDF

`predictive_pdf()` evaluates the predictive PDF of the model specified
by a
[`DataDistribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
`dist` and
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md) at
the given stage-one outcome.

## Usage

``` r
predictive_pdf(dist, prior, x1, n1, ...)

# S4 method for class 'DataDistribution,ContinuousPrior,numeric'
predictive_pdf(
  dist,
  prior,
  x1,
  n1,
  k = 10 * (prior@support[2] - prior@support[1]) + 1,
  ...
)

# S4 method for class 'DataDistribution,PointMassPrior,numeric'
predictive_pdf(dist, prior, x1, n1, ...)
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

`numeric`, value of the predictive PDF

## Examples

``` r
tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
predictive_pdf(Normal(), tmp, 2, 20)
#> [1] 0.2302199

predictive_pdf(Normal(), PointMassPrior(.3, 1), 1.5, 20) # ~.343
#> [1] 0.3426953
```
