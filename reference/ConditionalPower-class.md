# (Conditional) Power of a Design

This score evaluates P\[X₂ \> c2(design, X₁) \| X₁ = x₁\]. Note that the
distribution of X₂ is the posterior predictive after observing X₁ = x₁.

## Usage

``` r
ConditionalPower(dist, prior, label = "Pr[x2>=c2(x1)|x1]")

Power(dist, prior, label = "Pr[x2>=c2(x1)]")

# S4 method for class 'ConditionalPower,TwoStageDesign'
evaluate(s, design, x1, optimization = FALSE, ...)
```

## Arguments

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- prior:

  a [`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)
  object

- label:

  object label (string)

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

[`Scores`](https://optad.github.io/adoptr/reference/Scores.md)

## Examples

``` r
prior <- PointMassPrior(.4, 1)
cp <- ConditionalPower(Normal(), prior)
evaluate(
   cp,
   TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L),
   x1 = 1
)
#> [1] 0.5
# these two are equivalent:
expected(cp, Normal(), prior)
#> E[Pr[x2>=c2(x1)|x1]]<Normal<two-armed>;PointMass<0.40>> 
Power(Normal(), prior)
#> Pr[x2>=c2(x1)] 
```
