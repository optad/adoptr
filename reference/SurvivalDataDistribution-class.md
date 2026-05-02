# Log-rank test

Implements the normal approximation of the log-rank test statistic.

## Usage

``` r
Survival(event_rate, two_armed = TRUE)

# S4 method for class 'Survival'
quantile(x, probs, n, theta, ...)

# S4 method for class 'Survival,numeric'
simulate(object, nsim, n, theta, seed = NULL, ...)
```

## Arguments

- event_rate:

  probability that a subject will eventually have an event

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

  object of class `Survival`

- nsim:

  number of simulation runs

- seed:

  random seed

## Slots

- `event_rate`:

  cf. parameter 'event_rate'

## See also

see
[`probability_density_function`](https://optad.github.io/adoptr/reference/probability_density_function.md)
and
[`cumulative_distribution_function`](https://optad.github.io/adoptr/reference/cumulative_distribution_function.md)
to evaluate the pdf and the cdf, respectively.

## Examples

``` r
datadist <- Survival(event_rate=0.6, two_armed=TRUE)
```
