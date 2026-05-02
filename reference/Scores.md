# Scores

In `adoptr` scores are used to assess the performance of a design. This
can be done either conditionally on the observed stage-one outcome or
unconditionally. Consequently, score objects are either of class
`ConditionalScore` or `UnconditionalScore`.

## Usage

``` r
expected(s, data_distribution, prior, ...)

# S4 method for class 'ConditionalScore'
expected(s, data_distribution, prior, label = NA_character_, ...)

evaluate(s, design, ...)

# S4 method for class 'IntegralScore,TwoStageDesign'
evaluate(s, design, optimization = FALSE, subdivisions = 10000L, ...)
```

## Arguments

- s:

  `Score` object

- data_distribution:

  [`DataDistribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- prior:

  a [`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)
  object

- ...:

  further optional arguments

- label:

  object label (string)

- design:

  object

- optimization:

  logical, if `TRUE` uses a relaxation to real parameters of the
  underlying design; used for smooth optimization.

- subdivisions:

  maximal number of subdivisions when evaluating an integral score using
  adaptive quadrature (optimization = FALSE)

## Value

No return value. Generic description of class `Score`.

## Details

All scores can be evaluated on a design using the `evaluate` method.
Note that `evaluate` requires a third argument `x1` for conditional
scores (observed stage-one outcome). Any `ConditionalScore` can be
converted to a `UnconditionalScore` by forming its expected value using
`expected`. The returned unconditional score is of class
`IntegralScore`.

## See also

[`ConditionalPower`](https://optad.github.io/adoptr/reference/ConditionalPower-class.md),
[`ConditionalSampleSize`](https://optad.github.io/adoptr/reference/ConditionalSampleSize-class.md),
[`composite`](https://optad.github.io/adoptr/reference/composite.md)

## Examples

``` r
design <- TwoStageDesign(
  n1    = 25,
  c1f   = 0,
  c1e   = 2.5,
  n2    = 50,
  c2    = 1.96,
  order = 7L
)
prior <- PointMassPrior(.3, 1)

# conditional
cp <- ConditionalPower(Normal(), prior)
expected(cp, Normal(), prior)
#> E[Pr[x2>=c2(x1)|x1]]<Normal<two-armed>;PointMass<0.30>> 
evaluate(cp, design, x1 = .5)
#> [1] 0.3227581

# unconditional
power <- Power(Normal(), prior)
evaluate(power, design)
#> [1] 0.3269562
evaluate(power, design, optimization = TRUE) # use non-adaptive quadrature
#> [1] 0.3269562

```
