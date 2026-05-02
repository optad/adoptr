# Continuous univariate prior distributions

`ContinuousPrior` is a sub-class of
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)
implementing a generic representation of continuous prior distributions
over a compact interval on the real line.

## Usage

``` r
ContinuousPrior(
  pdf,
  support,
  order = 10,
  label = NA_character_,
  tighten_support = FALSE,
  check_normalization = TRUE
)
```

## Arguments

- pdf:

  vectorized univariate PDF function

- support:

  numeric vector of length two with the bounds of the compact interval
  on which the pdf is positive.

- order:

  `integer`, integration order of the employed Gaussian quadrature
  integration rule to evaluate scores. Automatically set to
  `length(n2_pivots)` if  
  `length(n2_pivots) == length(c2_pivots) > 1`, otherwise c2 and n2 are
  taken to be constant in stage-two and replicated to match the number
  of pivots specified by `order`

- label:

  object label (string)

- tighten_support:

  logical indicating if the support should be tightened

- check_normalization:

  logical indicating if it should be checked that `pdf` defines a
  density.

## Slots

- `pdf`:

  cf. parameter 'pdf'

- `support`:

  cf. parameter 'support'

- `pivots`:

  normalized pivots for integration rule (in \[-1, 1\]) the actual
  pivots are scaled to the support of the prior

- `weights`:

  weights of of integration rule at `pivots` for approximating integrals
  over `delta`

## See also

Discrete priors are supported via
[`PointMassPrior`](https://optad.github.io/adoptr/reference/PointMassPrior-class.md)

## Examples

``` r
ContinuousPrior(function(x) 2*x, c(0, 1))
#> ContinuousPrior<[0,1]> 
```
