# Create a collection of constraints

`subject_to(...)` can be used to generate an object of class
`ConstraintsCollection` from an arbitrary number of (un)conditional
constraints.

## Usage

``` r
subject_to(...)

# S4 method for class 'ConstraintsCollection,TwoStageDesign'
evaluate(s, design, optimization = FALSE, ...)
```

## Arguments

- ...:

  either constraint objects (for `subject_to` or optional arguments
  passed to `evaluate`)

- s:

  object of class `ConstraintCollection`

- design:

  object

- optimization:

  logical, if `TRUE` uses a relaxation to real parameters of the
  underlying design; used for smooth optimization.

## Value

an object of class `ConstraintsCollection`

## See also

`subject_to` is intended to be used for constraint specification the
constraints in
[`minimize`](https://optad.github.io/adoptr/reference/minimize.md).

## Examples

``` r
# define type one error rate and power
toer  <- Power(Normal(), PointMassPrior(0.0, 1))
power <- Power(Normal(), PointMassPrior(0.4, 1))

# create constrain collection
subject_to(
  toer  <= 0.025,
  power >= 0.9
)
#> An object of class "ConstraintsCollection"
#> Slot "unconditional_constraints":
#> [[1]]
#> Pr[x2>=c2(x1)] <= 0.025 
#> 
#> [[2]]
#> -Pr[x2>=c2(x1)]  <= -0.9 
#> 
#> 
#> Slot "conditional_constraints":
#> list()
#> 
```
