# Student's t data distribution

Implements exact t-distributions instead of a normal approximation

## Usage

``` r
Student(two_armed = TRUE)

# S4 method for class 'Student'
quantile(x, probs, n, theta, ...)

# S4 method for class 'Student,numeric'
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

  object of class `Student`

- nsim:

  number of simulation runs

- seed:

  random seed

## See also

see
[`probability_density_function`](https://optad.github.io/adoptr/reference/probability_density_function.md)
and
[`cumulative_distribution_function`](https://optad.github.io/adoptr/reference/cumulative_distribution_function.md)
to evaluate the pdf and the cdf, respectively.

## Examples

``` r
datadist <- Student(two_armed = TRUE)
```
