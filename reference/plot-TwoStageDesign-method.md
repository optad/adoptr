# Plot `TwoStageDesign` with optional set of conditional scores

This method allows to plot the stage-two sample size and decision
boundary functions of a chosen design.

## Usage

``` r
# S4 method for class 'TwoStageDesign'
plot(x, y = NULL, ..., rounded = TRUE, k = 100)
```

## Arguments

- x:

  design to plot

- y:

  not used

- ...:

  further named `ConditinonalScores` to plot for the design and/or
  further graphic parameters

- rounded:

  should n-values be rounded?

- k:

  number of points to use for plotting

## Value

a plot of the two-stage design

## Details

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)
and user-defined elements of the class
[`ConditionalScore`](https://optad.github.io/adoptr/reference/Scores.md).

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)

## Examples

``` r
design <- TwoStageDesign(50, 0, 2, 50, 2, 5)
cp     <- ConditionalPower(dist = Normal(), prior = PointMassPrior(.4, 1))
plot(design, "Conditional Power" = cp, cex.axis = 2)

```
