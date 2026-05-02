# Maximum Sample Size of a Design

This score evaluates `max(n(d))` for a design `d`.

## Usage

``` r
MaximumSampleSize(label = "max(n(x1))")

# S4 method for class 'MaximumSampleSize,TwoStageDesign'
evaluate(s, design, optimization = FALSE, ...)
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

- ...:

  further optional arguments

## See also

[Scores](https://optad.github.io/adoptr/reference/Scores.md) for general
scores and
[ConditionalSampleSize](https://optad.github.io/adoptr/reference/ConditionalSampleSize-class.md)
for evaluating the sample size point-wise.

## Examples

``` r
design <- TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L)
mss    <- MaximumSampleSize()
evaluate(mss, design)
#> [1] 100
```
