# (Conditional) Sample Size of a Design

This score simply evaluates `n(d, x1)` for a design `d` and the
first-stage outcome `x1`. The data distribution and prior are only
relevant when it is integrated.

## Usage

``` r
ConditionalSampleSize(label = "n(x1)")

ExpectedSampleSize(dist, prior, label = "E[n(x1)]")

ExpectedNumberOfEvents(dist, prior, label = "E[n(x1)]")

# S4 method for class 'ConditionalSampleSize,TwoStageDesign'
evaluate(s, design, x1, optimization = FALSE, ...)
```

## Arguments

- label:

  object label (string)

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- prior:

  a [`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)
  object

- s:

  [`Score`](https://optad.github.io/adoptr/reference/Scores.md) object

- design:

  object

- x1:

  stage-one test statistic

- optimization:

  logical, if `TRUE` uses a relaxation to real parameters of the
  underlying design; used for smooth optimization.

- ...:

  further optional arguments

## See also

[Scores](https://optad.github.io/adoptr/reference/Scores.md)

## Examples

``` r
design <- TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L)
prior  <- PointMassPrior(.4, 1)

css   <- ConditionalSampleSize()
evaluate(css, design, c(0, .5, 3))
#> [1] 100 100  50

ess   <- ExpectedSampleSize(Normal(), prior)
ene <- ExpectedNumberOfEvents(Survival(0.7), PointMassPrior(1.7, 1))

# those two are equivalent
evaluate(ess, design)
#> [1] 73.86249
evaluate(expected(css, Normal(), prior), design)
#> [1] 73.86249
```
