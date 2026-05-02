# Data distributions

`DataDistribution` is an abstract class used to represent the
distribution of a sufficient statistic `x` given a sample size `n` and a
single parameter value `theta`.

## Arguments

- x:

  outcome

- n:

  sample size

- theta:

  distribution parameter

- ...:

  further optional arguments

## Details

This abstraction layer allows the representation of t-distributions
(unknown variance), normal distribution (known variance), and normal
approximation of a binary endpoint. Currently, the two implemented
versions are
[`Normal-class`](https://optad.github.io/adoptr/reference/NormalDataDistribution-class.md)
and
[`Binomial-class`](https://optad.github.io/adoptr/reference/BinomialDataDistribution-class.md).

The logical option `two_armed` allows to decide whether a one-arm or a
two-arm (the default) design should be computed. In the case of a
two-arm design all sample sizes are per group.

## Slots

- `two_armed`:

  Logical that indicates if a two-arm design is assumed.

## Examples

``` r
normaldist   <- Normal(two_armed = FALSE)
binomialdist <- Binomial(rate_control = .25, two_armed = TRUE)
```
