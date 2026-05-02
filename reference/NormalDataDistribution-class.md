# Normal data distribution

Implements a normal data distribution for z-values given an observed
z-value and stage size. Standard deviation is 1 and mean θ √n where θ is
the standardized effect size. The option `two_armed` can be set to
decide whether a one-arm or a two-arm design should be computed.

## Usage

``` r
Normal(two_armed = TRUE)

# S4 method for class 'Normal'
quantile(x, probs, n, theta, ...)

# S4 method for class 'Normal,numeric'
simulate(object, nsim, n, theta, seed = NULL, ...)
```

## Arguments

- two_armed:

  logical indicating if a two-armed trial is regarded

- x:

  outcome

- probs:

  vector of probabilities

- n:

  sample size

- theta:

  distribution parameter

- ...:

  further optional arguments

- object:

  object of class `Normal`

- nsim:

  number of simulation runs

- seed:

  random seed

## Details

See
[`DataDistribution-class`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
for more details.

## See also

see
[`probability_density_function`](https://optad.github.io/adoptr/reference/probability_density_function.md)
and
[`cumulative_distribution_function`](https://optad.github.io/adoptr/reference/cumulative_distribution_function.md)
to evaluate the pdf and the cdf, respectively.

## Examples

``` r
datadist <- Normal(two_armed = TRUE)
```
