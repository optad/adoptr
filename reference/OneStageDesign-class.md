# One-stage designs

`OneStageDesign` implements a one-stage design as special case of a
two-stage design, i.e. as sub-class of
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md).
This is possible by defining n₂ = 0, c = c₁^(f) = c₁^(e), c₂(x₁) =
ifelse(x₁ \< c, Inf, -Inf). No integration pivots etc are required (set
to `NaN`).

## Usage

``` r
OneStageDesign(n, ...)

# S4 method for class 'numeric'
OneStageDesign(n, c, event_rate)

# S4 method for class 'OneStageDesign'
TwoStageDesign(n1, event_rate, order = 5L, eps = 0.01, ...)

# S4 method for class 'OneStageDesignSurvival'
TwoStageDesign(n1, order = 5L, eps = 0.01, ...)

# S4 method for class 'OneStageDesign'
plot(x, y, ...)
```

## Arguments

- n:

  sample size (stage-one sample size)

- ...:

  further optional arguments

- c:

  rejection boundary (c = c₁^(f) = c₁^(e))

- event_rate:

  probability that a subject in either group will eventually have an
  event, only needs to be specified for time-to-event endpoints.

- n1:

  `OneStageDesign` object to convert, overloaded from
  [`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)

- order:

  integer \>= 2, default is 5; order of Gaussian quadrature integration
  rule to use for new TwoStageDesign.

- eps:

  numeric \> 0, default = .01; the single critical value c must be split
  in a continuation interval \[c1f, c1e\]; this is given by c +/- eps.

- x:

  design to plot

- y:

  not used

## Details

Note that the default
[`plot,TwoStageDesign-method`](https://optad.github.io/adoptr/reference/plot-TwoStageDesign-method.md)
method is not supported for `OneStageDesign` objects.

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md),
[`GroupSequentialDesign-class`](https://optad.github.io/adoptr/reference/GroupSequentialDesign-class.md)

## Examples

``` r
design <- OneStageDesign(30, 1.96)
summary(design)
#> OneStageDesign: n1 =  30 
#>           futility | continue | efficacy
#>       x1:     1.96 |   NaN |  1.96
#>   c2(x1):     +Inf |    NA |  -Inf
#>   n2(x1):        0 |     0 |     0
#> 
design_twostage <- TwoStageDesign(design)
summary(design_twostage)
#> TwoStageDesign: n1 =  30 
#>           futility |            continue           | efficacy
#>       x1:     1.95 |  1.95  1.95  1.96  1.97  1.97 |  1.97
#>   c2(x1):     +Inf | +3.00 +3.00 +0.00 -3.00 -3.00 |  -Inf
#>   n2(x1):        0 |     0     0     0     0     0 |     0
#> 
design_survival <- OneStageDesign(30, 1.96, 0.7)

TwoStageDesign(design_survival)
#> TwoStageDesignSurvival<n_events1=30;1.9<=x1<=2.0;n_events2=0> 
```
