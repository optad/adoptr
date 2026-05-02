# Fix parameters during optimization

The methods `make_fixed` and `make_tunable` can be used to modify the
'tunability' status of parameters in a
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)
object. Tunable parameters are optimized over, non-tunable ('fixed')
parameters are considered given and not altered during optimization.

## Usage

``` r
make_tunable(x, ...)

# S4 method for class 'TwoStageDesign'
make_tunable(x, ...)

make_fixed(x, ...)

# S4 method for class 'TwoStageDesign'
make_fixed(x, ...)
```

## Arguments

- x:

  `TwoStageDesign` object

- ...:

  unquoted names of slots for which the tunability status should be
  changed.

## Value

an updated object of class
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md),
[`tunable_parameters`](https://optad.github.io/adoptr/reference/tunable_parameters.md)
for converting tunable parameters of a design object to a numeric vector
(and back), and
[`minimize`](https://optad.github.io/adoptr/reference/minimize.md) for
the actual minimzation procedure

## Examples

``` r
design <- TwoStageDesign(25, 0, 2, 25, 2, order = 5)
# default: all parameters are tunable (except integration pivots,
# weights and tunability status itself)
design@tunable
#>             n1            c1f            c1e      n2_pivots      c2_pivots 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> x1_norm_pivots        weights        tunable 
#>          FALSE          FALSE          FALSE 

# make n1 and the pivots of n2 fixed (not changed during optimization)
design <- make_fixed(design, n1, n2_pivots)
design@tunable
#>             n1            c1f            c1e      n2_pivots      c2_pivots 
#>          FALSE           TRUE           TRUE          FALSE           TRUE 
#> x1_norm_pivots        weights        tunable 
#>          FALSE          FALSE          FALSE 

# make them tunable again
design <- make_tunable(design, n1, n2_pivots)
design@tunable
#>             n1            c1f            c1e      n2_pivots      c2_pivots 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> x1_norm_pivots        weights        tunable 
#>          FALSE          FALSE          FALSE 
```
