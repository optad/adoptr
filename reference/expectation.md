# Expected value of a function

Computes the expected value of a vectorized, univariate function `f`
with respect to a distribution `dist`. I.e., E\[f(X)\].

## Usage

``` r
expectation(dist, f, ...)

# S4 method for class 'ContinuousPrior,function'
expectation(dist, f, ...)

# S4 method for class 'PointMassPrior,function'
expectation(dist, f, ...)
```

## Arguments

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- f:

  a univariate function, must be vectorized

- ...:

  further optional arguments

## Value

`numeric`, expected value of `f` with respect to `dist`

## Examples

``` r
expectation(
    ContinuousPrior(function(x) dunif(x, .2, .4), c(.2, .4)),
    identity
)
#> [1] 0.3
# > 0.3

expectation(PointMassPrior(c(0, .5), c(.3, .7)), identity)
#> [1] 0.35
# > .35
```
