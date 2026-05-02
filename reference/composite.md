# Score Composition

`composite` defines new composite scores by point-wise evaluation of
scores in any valid numerical expression.

## Usage

``` r
composite(expr, label = NA_character_)

# S4 method for class 'CompositeScore,TwoStageDesign'
evaluate(s, design, ...)
```

## Arguments

- expr:

  Expression (in curly brackets); must contain at least one score
  variable; if multiple scores are used, they must either all be
  conditional or unconditional. Currently, no non-score variables are
  supported

- label:

  object label (string)

- s:

  object of class `CompositeScore`

- design:

  object

- ...:

  further optional arguments

## Value

an object of class `CompositeConditionalScore` or
`CompositeUnconditionalScore` depending on the class of the scores used
in `expr`

## See also

[Scores](https://optad.github.io/adoptr/reference/Scores.md)

## Examples

``` r
ess   <- ExpectedSampleSize(Normal(), PointMassPrior(.4, 1))
power <- Power(Normal(), PointMassPrior(.4, 1))

# linear combination:
composite({ess - 50*power})
#> E[n(x1)]  - 50 * Pr[x2>=c2(x1)]  

# control flow (e.g. for and while loops)
composite({
  res <- 0
  for (i in 1:3) {
     res <- res + ess
  }
  res
})
#> res <- 0; for (i in 1:3) {
#>     res <- res + E[n(x1)] 
#> }; res 

# functional composition
composite({log(ess)})
#> log(E[n(x1)] ) 
cp <- ConditionalPower(Normal(), PointMassPrior(.4, 1))
composite({3*cp})
#> 3 * Pr[x2>=c2(x1)|x1]  
```
