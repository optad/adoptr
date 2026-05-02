# Composite Scores

``` r

library(adoptr)
```

While adopt also allows implementation of custom scores via subclassing,
for most applications a simple point-wise arithmetic on scores is
sufficient. For instance, consider the case of a utility maximizing
approach to planning where not a hard constraint on power but rather a
trade-off between power and expected sample size is required. The
simplest utility function would just be a weighted sum of both power
(negative weight since we minimize costs!) and expected sample size.

Consider the following situation

``` r

H_0      <- PointMassPrior(.0, 1)
H_1      <- PointMassPrior(.2, 1)
datadist <- Binomial(.1, two_armed = FALSE)

ess   <- ExpectedSampleSize(datadist, H_1)
power <- Power(datadist, H_1)
toer  <- Power(datadist, H_0)
```

Adoptr supports such `CompositeScores` via the `composite` function:

``` r

objective <- composite({ess - 50*power})
```

The new unconditional score can be evaluated as usual, e.g.

``` r

design <- TwoStageDesign(
    n1  = 100,
    c1f = .0,
    c1e = 2.0,
    n2_pivots = rep(150, 5),
    c2_pivots = sapply(1 + adoptr:::GaussLegendreRule(5)$nodes, function(x) -x + 2)
)

evaluate(objective, design)
#> [1] 50.16813
```

Note that conditional and unconditional scores cannot be mixed in an
expression passed to `composite`. Composite conditional score, however,
are possible as well.

``` r

cp  <- ConditionalPower(datadist, H_1)
css <- ConditionalSampleSize()

cs  <- composite({css - 50*cp})
```

``` r

evaluate(cs, design, c(0, .5, 1))
#> [1] 200.0014 200.0003 200.0001
```

Of course, composite conditional scores can also be integrated

``` r

evaluate(expected(cs, datadist, H_1), design)
#> [1] 50.16813
```

and (due to linearity) the result is exactly the same as before.

## Functional Composition

Composite scores are not restricted to linear operations but support any
valid numerical expression:

``` r

cs <- composite({log(css) - 50*sin(cp)})
evaluate(cs, design, c(0, .5, 1))
#> [1] -36.55135 -36.55192 -36.55205
```

Even control flow is supported:

``` r

cs <- composite({
  res <- 0
  for (i in 1:3) {
    res <- res + css
  }
  res
})
evaluate(cs, design, c(0, .5, 1))
#> [1] 750 750 750
```

The only real constraint is that the expression must be vectorized.
