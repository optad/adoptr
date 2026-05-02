# Condition a prior on an interval

Restrict an object of class
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md) to a
sub-interval and re-normalize the PDF.

## Usage

``` r
condition(dist, interval, ...)

# S4 method for class 'ContinuousPrior,numeric'
condition(dist, interval, ...)

# S4 method for class 'PointMassPrior,numeric'
condition(dist, interval, ...)
```

## Arguments

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- interval:

  length-two numeric vector giving the parameter interval to condition
  on

- ...:

  further optional arguments

## Value

conditional
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md) on
given interval

## Examples

``` r
tmp <- condition(
    ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4)),
    c(.3, .5)
)
bounds(tmp) # c(.3, .4)
#> [1] 0.3 0.4

tmp <- condition(PointMassPrior(c(0, .5), c(.3, .7)), c(-1, .25))
expectation(tmp, identity) # 0
#> [1] 0
```
