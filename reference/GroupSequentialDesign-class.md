# Group-sequential two-stage designs

Group-sequential designs are a sub-class of the `TwoStageDesign` class
with constant stage-two sample size. See
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)
for slot details. Any group-sequential design can be converted to a
fully flexible `TwoStageDesign` (see examples section).

## Usage

``` r
GroupSequentialDesign(n1, ...)

# S4 method for class 'numeric'
GroupSequentialDesign(
  n1,
  c1f,
  c1e,
  n2_pivots,
  c2_pivots,
  order = NULL,
  event_rate,
  ...
)

# S4 method for class 'GroupSequentialDesign'
TwoStageDesign(n1, event_rate, ...)

# S4 method for class 'GroupSequentialDesignSurvival'
TwoStageDesign(n1, ...)
```

## Arguments

- n1:

  stage one sample size or `GroupSequentialDesign` object to convert
  (overloaded from
  [`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md))

- ...:

  further optional arguments

- c1f:

  early futility stopping boundary

- c1e:

  early efficacy stopping boundary

- n2_pivots:

  numeric of length one, stage-two sample size

- c2_pivots:

  numeric vector, stage-two critical values on the integration pivot
  points

- order:

  of the Gaussian quadrature rule to use for integration, set to
  length(c2_pivots) if NULL, otherwise first value of c2_pivots is
  repeated 'order'-times.

- event_rate:

  probability that a subject in either group will eventually have an
  event, only needs to be specified for time-to-event endpoints.

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)
for superclass and inherited methods

## Examples

``` r
design <- GroupSequentialDesign(25, 0, 2, 25, c(1, 1.5, 2.5))
summary(design)
#> GroupSequentialDesign: n1 =  25 
#>           futility |      continue     | efficacy
#>       x1:    -0.00 |  0.23  1.00  1.77 |  2.00
#>   c2(x1):     +Inf | +1.00 +1.50 +2.50 |  -Inf
#>   n2(x1):        0 |    25    25    25 |     0
#> 

design_survival <- GroupSequentialDesign(25, 0, 2, 25, c(1, 1.5, 2.5), event_rate = 0.7)

TwoStageDesign(design)
#> TwoStageDesign<n1=25;0.0<=x1<=2.0:n2=25> 

TwoStageDesign(design_survival)
#> TwoStageDesignSurvival<n_events1=25;0.0<=x1<=2.0;n_events2=25> 
```
