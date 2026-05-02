# Get support of a prior or data distribution

`bounds()` returns the range of the support of a prior or data
distribution.

## Usage

``` r
bounds(dist, ...)

# S4 method for class 'ContinuousPrior'
bounds(dist, ...)

# S4 method for class 'PointMassPrior'
bounds(dist, ...)
```

## Arguments

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- ...:

  further optional arguments

## Value

`numeric` of length two, `c(lower, upper)`

## Examples

``` r
bounds(ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4)))
#> [1] 0.2 0.4
# > 0.2 0.4

bounds(PointMassPrior(c(0, .5), c(.3, .7)))
#> [1] 0.0 0.5
# > 0.3 0.7
```
