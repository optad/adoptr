# Query sample size of a design

Methods to access the stage-one, stage-two, or overall sample size of a
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md).
`n1` returns the first-stage sample size of a design, `n2` the stage-two
sample size conditional on the stage-one test statistic and `n` the
overall sample size `n1 + n2`. Internally, objects of the class
`TwoStageDesign` allow non-natural, real sample sizes to allow smooth
optimization (cf.
[`minimize`](https://optad.github.io/adoptr/reference/minimize.md) for
details). The optional argument `round` allows to switch between the
internal real representation and a rounded version (rounding to the next
positive integer).

## Usage

``` r
n1(d, ...)

# S4 method for class 'TwoStageDesign'
n1(d, round = TRUE, ...)

n2(d, x1, ...)

# S4 method for class 'TwoStageDesign,numeric'
n2(d, x1, round = TRUE, ...)

n(d, x1, ...)

# S4 method for class 'TwoStageDesign,numeric'
n(d, x1, round = TRUE, ...)

# S4 method for class 'OneStageDesign,numeric'
n2(d, x1, ...)

# S4 method for class 'GroupSequentialDesign,numeric'
n2(d, x1, round = TRUE, ...)
```

## Arguments

- d:

  design

- ...:

  further optional arguments

- round:

  `logical` should sample sizes be rounded to next integer?

- x1:

  stage-one test statistic

## Value

sample size value of design `d` at point `x1`

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md),
see [`c2`](https://optad.github.io/adoptr/reference/critical-values.md)
for accessing the critical values

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

n1(design) # 25
#> [1] 25
design@n1 # 25
#> [1] 25

n(design, x1 = 2.2) # 75
#> [1] 75

```
