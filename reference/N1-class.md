# Regularize n1

`N1` is a class that computes the `n1` value of a design. This can be
used as a score in
[`minimize`](https://optad.github.io/adoptr/reference/minimize.md).

## Usage

``` r
N1(label = NA_character_)

# S4 method for class 'N1,TwoStageDesign'
evaluate(s, design, optimization = FALSE, ...)
```

## Arguments

- label:

  object label (string)

- s:

  [`Score`](https://optad.github.io/adoptr/reference/Scores.md) object

- design:

  object

- optimization:

  logical, if `TRUE` uses a relaxation to real parameters of the
  underlying design; used for smooth optimization.

- ...:

  further optional arguments

## Value

an object of class `N1`

## See also

See
[`AverageN2`](https://optad.github.io/adoptr/reference/AverageN2-class.md)
for a regularization of the second-stage sample size.

## Examples

``` r
n1_score <- N1()

evaluate(
   N1(),
   TwoStageDesign(70, 0, 2, rep(60, 6), rep(1.7, 6))
) # 70
#> [1] 70
```
