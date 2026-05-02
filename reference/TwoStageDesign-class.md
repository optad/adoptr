# Two-stage designs

`TwoStageDesign` is the fundamental design class of the
[adoptr](https://optad.github.io/adoptr/reference/adoptr.md) package.
Formally, we represent a generic two-stage design as a five-tuple (n₁,
c₁^(f), c₁^(e), n₂(·), c₂(·)). Here, n₁ is the first-stage sample size
(per group), c₁^(f) and c₁^(e) are boundaries for early stopping for
futility and efficacy, respectively. Since the trial design is a
two-stage design, the elements n₂(·) (stage-two sample size) and c₂(·)
(stage-two critical value) are functions of the first-stage outcome
X₁=x₁. X₁ denotes the first-stage test statistic. A brief description on
this definition of two-stage designs can be read
[here](https://optad.github.io/adoptr/articles/adoptr.html). For
available methods, see the 'See Also' section at the end of this page.

## Usage

``` r
TwoStageDesign(n1, ...)

# S4 method for class 'numeric'
TwoStageDesign(
  n1,
  c1f,
  c1e,
  n2_pivots,
  c2_pivots,
  order = NULL,
  event_rate,
  ...
)

# S4 method for class 'TwoStageDesign'
summary(object, ..., rounded = TRUE)
```

## Arguments

- n1:

  stage-one sample size

- ...:

  further optional arguments

- c1f:

  early futility stopping boundary

- c1e:

  early efficacy stopping boundary

- n2_pivots:

  numeric vector, stage-two sample size on the integration pivot points

- c2_pivots:

  numeric vector, stage-two critical values on the integration pivot
  points

- order:

  `integer`, integration order of the employed Gaussian quadrature
  integration rule to evaluate scores. Automatically set to
  `length(n2_pivots)` if  
  `length(n2_pivots) == length(c2_pivots) > 1`, otherwise c2 and n2 are
  taken to be constant in stage-two and replicated to match the number
  of pivots specified by `order`

- event_rate:

  probability that a subject in either group will eventually have an
  event, only needs to be specified for time-to-event endpoints

- object:

  object to show

- rounded:

  should rounded n-values be used?

## Details

`summary` can be used to quickly compute and display basic facts about a
TwoStageDesign. An arbitrary number of names
[`UnconditionalScore`](https://optad.github.io/adoptr/reference/Scores.md)
objects can be provided via the optional arguments `...` and are
included in the summary displayed using
[`print`](https://optad.github.io/adoptr/reference/print.adoptrOptimizationResult.md).

## Slots

- `n1`:

  cf. parameter 'n1'

- `c1f`:

  cf. parameter 'c1f'

- `c1e`:

  cf. parameter 'c1e'

- `n2_pivots`:

  vector of length 'order' giving the values of n2 at the pivot points
  of the numeric integration rule

- `c2_pivots`:

  vector of length order giving the values of c2 at the pivot points of
  the numeric integration rule

- `x1_norm_pivots`:

  normalized pivots for integration rule (in \[-1, 1\]) the actual
  pivots are scaled to the interval \[c1f, c1e\] and can be obtained by
  the internal method  
  `adoptr:::scaled_integration_pivots(design)`

- `weights`:

  weights of of integration rule at `x1_norm_pivots` for approximating
  integrals over `x1`

- `tunable`:

  named logical vector indicating whether corresponding slot is
  considered a tunable parameter (i.e. whether it can be changed during
  optimization via
  [`minimize`](https://optad.github.io/adoptr/reference/minimize.md) or
  not; cf.  
  [`make_fixed`](https://optad.github.io/adoptr/reference/make_tunable.md))

## See also

For accessing sample sizes and critical values safely, see methods in
[`n`](https://optad.github.io/adoptr/reference/n.md) and
[`c2`](https://optad.github.io/adoptr/reference/critical-values.md); for
modifying behaviour during optimizaton see
[`make_tunable`](https://optad.github.io/adoptr/reference/make_tunable.md);
to convert between S4 class represenation and numeric vector, see
[`tunable_parameters`](https://optad.github.io/adoptr/reference/tunable_parameters.md);
for simulating from a given design, see
[`simulate`](https://optad.github.io/adoptr/reference/simulate-TwoStageDesign-numeric-method.md);
for plotting see
[`plot,TwoStageDesign-method`](https://optad.github.io/adoptr/reference/plot-TwoStageDesign-method.md).
Both
[group-sequential](https://optad.github.io/adoptr/reference/GroupSequentialDesign-class.md)
and [one-stage
designs](https://optad.github.io/adoptr/reference/OneStageDesign-class.md)
(!) are implemented as subclasses of `TwoStageDesign`.

## Examples

``` r
design <- TwoStageDesign(50, 0, 2, 50.0, 2.0, 5)
pow    <- Power(Normal(), PointMassPrior(.4, 1))
summary(design, "Power" = pow)
#> TwoStageDesign: n1 =  50 
#>           futility |            continue           | efficacy
#>       x1:    -0.00 |  0.09  0.46  1.00  1.54  1.91 |  2.00
#>   c2(x1):     +Inf | +2.00 +2.00 +2.00 +2.00 +2.00 |  -Inf
#>   n2(x1):        0 |    50    50    50    50    50 |     0
#>    Power:      0.739
#> 
```
