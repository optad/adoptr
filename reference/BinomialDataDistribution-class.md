# Binomial data distribution

Implements the normal approximation for a test on rates. The reponse
rate in the control group, r_(C), has to be specified by `rate_control`.
The null hypothesis is: r_(E) ≤ r_(C), where r_(E) denotes the response
rate in the invervention group. It is tested against the alternative
r_(E) \> r_(C). The test statistic is given as X₁ = √n (r_(E) - r_(C)) /
√(2 r₀ (1-r₀)), where r₀ denotes the mean between r_(E) and r_(C) in the
two-armed case, and r_(E) in the one-armed case.#' All priors have to be
defined for the rate difference r_(E) - r_(C).

## Usage

``` r
Binomial(rate_control, two_armed = TRUE)

# S4 method for class 'Binomial'
quantile(x, probs, n, theta, ...)

# S4 method for class 'Binomial,numeric'
simulate(object, nsim, n, theta, seed = NULL, ...)
```

## Arguments

- rate_control:

  assumed response rate in control group

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

  object of class `Binomial`

- nsim:

  number of simulation runs

- seed:

  random seed

## Details

Note that `simulate` for class `Binomial` simulates the normal
approximation of the test statistic.

## Slots

- `rate_control`:

  cf. parameter 'rate_control'

## See also

see
[`probability_density_function`](https://optad.github.io/adoptr/reference/probability_density_function.md)
and
[`cumulative_distribution_function`](https://optad.github.io/adoptr/reference/cumulative_distribution_function.md)
to evaluate the pdf and the cdf, respectively.

## Examples

``` r
datadist <- Binomial(rate_control = 0.2, two_armed = FALSE)
```
