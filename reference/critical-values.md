# Query critical values of a design

Methods to access the stage-two critical values of a
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md).
`c2` returns the stage-two critical value conditional on the stage-one
test statistic.

## Usage

``` r
c2(d, x1, ...)

# S4 method for class 'TwoStageDesign,numeric'
c2(d, x1, ...)

# S4 method for class 'OneStageDesign,numeric'
c2(d, x1, ...)
```

## Arguments

- d:

  design

- x1:

  stage-one test statistic

- ...:

  further optional arguments

## Value

the critical value function `c2` of design `d` at position `x1`

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md),
see [`n`](https://optad.github.io/adoptr/reference/n.md) for accessing
the sample size of a design

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

c2(design, 2.2) # 1.96
#> [1] 1.96
c2(design, 3.0) # -Inf
#> [1] -Inf
c2(design, -1.0) # Inf
#> [1] Inf

design <- TwoStageDesign(
   n1    = 25,
   c1f   = 0,
   c1e   = 2.5,
   n2    = 50,
   c2    = 1.96,
   order = 7L
)

c2(design, 2.2) # 1.96
#> [1] 1.96
c2(design, 3.0) # -Inf
#> [1] -Inf
c2(design, -1.0) # Inf
#> [1] Inf
```
