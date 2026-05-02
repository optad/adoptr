# Switch between numeric and S4 class representation of a design

Get tunable parameters of a design as numeric vector via
`tunable_parameters` or `update` a design object with a suitable vector
of values for its tunable parameters.

## Usage

``` r
tunable_parameters(object, ...)

# S4 method for class 'TwoStageDesign'
tunable_parameters(object, ...)

# S4 method for class 'TwoStageDesign'
update(object, params, ...)

# S4 method for class 'OneStageDesign'
update(object, params, ...)
```

## Arguments

- object:

  `TwoStageDesign` object to update

- ...:

  further optional arguments

- params:

  vector of design parameters, must be in same order as returned by  
  `tunable_parameters`

## Value

`tunable_parameters` returns the numerical values of all tunable
parameters as a vector. `update` returns the updated design.

## Details

The `tunable` slot of a
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)
stores information about the set of design parameters which are
considered fixed (not changed during optimization) or tunable (changed
during optimization). For details on how to fix certain parameters or
how to make them tunable again, see
[`make_fixed`](https://optad.github.io/adoptr/reference/make_tunable.md)
and
[`make_tunable`](https://optad.github.io/adoptr/reference/make_tunable.md).

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)

## Examples

``` r
design  <- TwoStageDesign(25, 0, 2, 25, 2, order = 5)
tunable_parameters(design)
#>  [1] 25  0  2 25 25 25 25 25  2  2  2  2  2
design2 <- update(design, tunable_parameters(design) + 1)
tunable_parameters(design2)
#>  [1] 26  1  3 26 26 26 26 26  3  3  3  3  3
```
