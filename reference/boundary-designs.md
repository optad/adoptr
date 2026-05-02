# Boundary designs

The optimization method
[`minimize`](https://optad.github.io/adoptr/reference/minimize.md) is
based on the package `nloptr`. This requires upper and lower boundaries
for optimization. Such boundaries can be computed via
`lower_boundary_design` respectively `upper_boundary_design`. They are
implemented by default in
[`minimize`](https://optad.github.io/adoptr/reference/minimize.md). Note
that [`minimize`](https://optad.github.io/adoptr/reference/minimize.md)
allows the user to define its own boundary designs, too.

## Usage

``` r
get_lower_boundary_design(initial_design, ...)

get_upper_boundary_design(initial_design, ...)

# S4 method for class 'OneStageDesign'
get_lower_boundary_design(initial_design, n1 = 1, c1_buffer = 2, ...)

# S4 method for class 'GroupSequentialDesign'
get_lower_boundary_design(
  initial_design,
  n1 = 1,
  n2_pivots = 1,
  c1_buffer = 2,
  c2_buffer = 2,
  ...
)

# S4 method for class 'TwoStageDesign'
get_lower_boundary_design(
  initial_design,
  n1 = 1,
  n2_pivots = 1,
  c1_buffer = 2,
  c2_buffer = 2,
  ...
)

# S4 method for class 'OneStageDesign'
get_upper_boundary_design(
  initial_design,
  n1 = 5 * initial_design@n1,
  c1_buffer = 2,
  ...
)

# S4 method for class 'GroupSequentialDesign'
get_upper_boundary_design(
  initial_design,
  n1 = 5 * initial_design@n1,
  n2_pivots = 5 * initial_design@n2_pivots,
  c1_buffer = 2,
  c2_buffer = 2,
  ...
)

# S4 method for class 'TwoStageDesign'
get_upper_boundary_design(
  initial_design,
  n1 = 5 * initial_design@n1,
  n2_pivots = 5 * initial_design@n2_pivots,
  c1_buffer = 2,
  c2_buffer = 2,
  ...
)
```

## Arguments

- initial_design:

  The initial design

- ...:

  optional arguments

  The values `c1f` and `c1e` from the initial design are shifted to
  `c1f - c1_buffer` and `c1e - c1_buffer` in
  `get_lower_boundary_design`, respectively, to  
  `c1f + c1_buffer` and `c1e + c1_buffer` in
  `get_upper_boundary_design`. This is handled analogously with
  `c2_pivots` and `c2_buffer`.

- n1:

  bound for the first-stage sample size n1

- c1_buffer:

  shift of the early-stopping boundaries from the initial ones

- n2_pivots:

  bound for the second-stage sample size n2

- c2_buffer:

  shift of the final decision boundary from the initial one

## Value

An object of class
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md).

## Examples

``` r
initial_design <- TwoStageDesign(
  n1    = 25,
  c1f   = 0,
  c1e   = 2.5,
  n2    = 50,
  c2    = 1.96,
  order = 7L
  )
get_lower_boundary_design(initial_design)
#> TwoStageDesign<n1=1;-2.0<=x1<=0.5:n2=1> 
```
