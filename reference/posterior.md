# Compute posterior distribution

Return posterior distribution given observing stage-one outcome.

## Usage

``` r
posterior(dist, prior, x1, n1, ...)

# S4 method for class 'DataDistribution,ContinuousPrior,numeric'
posterior(dist, prior, x1, n1, ...)

# S4 method for class 'DataDistribution,PointMassPrior,numeric'
posterior(dist, prior, x1, n1, ...)
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

## Value

Object of class
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)

## Examples

``` r
tmp <- ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4))
posterior(Normal(), tmp, 2, 20)
#> ContinuousPrior<[0.2,0.4]> 

posterior(Normal(), PointMassPrior(0, 1), 2, 20)
#> PointMass<0.00> 
```
