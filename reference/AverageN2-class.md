# Regularization via L1 norm

Implements the L1-norm of the design's stage-two sample size function.
The average of the stage-two sample size without weighting with the data
distribution is computed. This can be interpreted as integration over a
unifrom prior on the continuation region.

## Usage

``` r
AverageN2(label = NA_character_)

# S4 method for class 'AverageN2,TwoStageDesign'
evaluate(s, design, optimization = FALSE, subdivisions = 10000L, ...)
```

## Arguments

- label:

  object label (string)

- s:

  [`Score`](https://optad.github.io/adoptr/reference/Scores.md) object

- design:

  object

- optimization:

  logical, if `TRUE` uses a relaxation to real parameters of the
  underlying design; used for smooth optimization.

- subdivisions:

  number of subdivisions to use for adaptive integration (only affects
  non-optimization code)

- ...:

  further optional arguments

## Value

an object of class `AverageN2`

## See also

[`N1`](https://optad.github.io/adoptr/reference/N1-class.md) for
penalizing n1 values

## Examples

``` r
avn2 <- AverageN2()

evaluate(
   AverageN2(),
   TwoStageDesign(100, 0.5, 1.5, 60.0, 1.96, order = 5L)
) # 60
#> [1] 60
```
