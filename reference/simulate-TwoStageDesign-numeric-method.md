# Draw samples from a two-stage design

`simulate` allows to draw samples from a given
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md).

## Usage

``` r
# S4 method for class 'TwoStageDesign,numeric'
simulate(object, nsim, dist, theta, seed = NULL, ...)
```

## Arguments

- object:

  `TwoStageDesign` to draw samples from

- nsim:

  number of simulation runs

- dist:

  data distribution

- theta:

  location parameter of the data distribution

- seed:

  random seed

- ...:

  further optional arguments

## Value

[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns a
`data.frame` with `nsim` rows and for each row (each simulation run) the
following columns

- theta: The effect size

- n1: First-stage sample size

- c1f: Stopping for futility boundary

- c1e: Stopping for efficacy boundary

- x1: First-stage outcome

- n2: Resulting second-stage sample size after observing x1

- c2: Resulting second-stage decision-boundary after observing x1

- x2: Second-stage outcome

- reject: Decision whether the null hypothesis is rejected or not

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)

## Examples

``` r
design <- TwoStageDesign(25, 0, 2, 25, 2, order = 5)
# draw samples assuming two-armed design
simulate(design, 10, Normal(), .3, 42)
#>    theta n1 c1f c1e        x1 n2   c2         x2 reject
#> 1    0.3 25   0   2 2.4316186  0 -Inf  1.3048697   TRUE
#> 2    0.3 25   0   2 0.4959620 25    2  3.3473056   TRUE
#> 3    0.3 25   0   2 1.4237886 25    2 -0.3282005  FALSE
#> 4    0.3 25   0   2 1.6935228 25    2  0.7818714  FALSE
#> 5    0.3 25   0   2 1.4649285 25    2  0.9273388  FALSE
#> 6    0.3 25   0   2 0.9545357 25    2  1.6966106  FALSE
#> 7    0.3 25   0   2 2.5721822  0 -Inf -0.2842529   TRUE
#> 8    0.3 25   0   2 0.9660011 25    2 -1.5957952  FALSE
#> 9    0.3 25   0   2 3.0790839  0 -Inf -2.4404669   TRUE
#> 10   0.3 25   0   2 0.9979461 25    2  2.3807735   TRUE
```
