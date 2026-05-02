# Univariate discrete point mass priors

`PointMassPrior` is a sub-class of
[`Prior`](https://optad.github.io/adoptr/reference/Prior-class.md)
representing a univariate prior over a discrete set of points with
positive probability mass.

## Usage

``` r
PointMassPrior(theta, mass, label = NA_character_)
```

## Arguments

- theta:

  numeric vector of pivot points with positive prior mass

- mass:

  numeric vector of probability masses at the pivot points (must sum to
  1)

- label:

  object label (string)

## Value

an object of class `PointMassPrior`, `theta` is automatically sorted in
ascending order

## Slots

- `theta`:

  cf. parameter 'theta'

- `mass`:

  cf. parameter 'mass'

## See also

To represent continuous prior distributions use
[`ContinuousPrior`](https://optad.github.io/adoptr/reference/ContinuousPrior-class.md).

## Examples

``` r
PointMassPrior(c(0, .5), c(.3, .7))
#> PointMass<Pr[0.00]=0.30;Pr[0.50]=0.70> 
```
