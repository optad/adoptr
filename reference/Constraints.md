# Formulating Constraints

Conceptually, constraints work very similar to scores (any score can be
put in a constraint). Currently, constraints of the form 'score \<=/\>=
x', 'x \<=/\>= score' and 'score \<=/\>= score' are admissible.

## Usage

``` r
# S4 method for class 'Constraint,TwoStageDesign'
evaluate(s, design, optimization = FALSE, ...)

# S4 method for class 'ConditionalScore,numeric'
e1 <= e2

# S4 method for class 'ConditionalScore,numeric'
e1 >= e2

# S4 method for class 'numeric,ConditionalScore'
e1 <= e2

# S4 method for class 'numeric,ConditionalScore'
e1 >= e2

# S4 method for class 'ConditionalScore,ConditionalScore'
e1 <= e2

# S4 method for class 'ConditionalScore,ConditionalScore'
e1 >= e2

# S4 method for class 'UnconditionalScore,numeric'
e1 <= e2

# S4 method for class 'UnconditionalScore,numeric'
e1 >= e2

# S4 method for class 'numeric,UnconditionalScore'
e1 <= e2

# S4 method for class 'numeric,UnconditionalScore'
e1 >= e2

# S4 method for class 'UnconditionalScore,UnconditionalScore'
e1 <= e2

# S4 method for class 'UnconditionalScore,UnconditionalScore'
e1 >= e2
```

## Arguments

- s:

  [`Score`](https://optad.github.io/adoptr/reference/Scores.md) object

- design:

  object

- optimization:

  logical, if `TRUE` uses a relaxation to real parameters of the
  underlying design; used for smooth optimization.

- ...:

  further optional arguments

- e1:

  left hand side (score or numeric)

- e2:

  right hand side (score or numeric)

## Value

an object of class `Constraint`

## See also

[`minimize`](https://optad.github.io/adoptr/reference/minimize.md)

## Examples

``` r
design <- OneStageDesign(50, 1.96)

cp     <- ConditionalPower(Normal(), PointMassPrior(0.4, 1))
pow    <- Power(Normal(), PointMassPrior(0.4, 1))

# unconditional power constraint
constraint1 <- pow >= 0.8
evaluate(constraint1, design)
#> [1] 0.2840466

# conditional power constraint
constraint2 <- cp  >= 0.7
evaluate(constraint2, design, .5)
#> [1] 0.7
constraint3 <- 0.7 <= cp # same as constraint2
evaluate(constraint3, design, .5)
#> [1] 0.7
```
