# Chi-Squared data distribution

Implements a chi-squared distribution. The classes `Pearson2xk` and
`ZSquared` are subclasses, used in two different situations.
`Pearson2xK` is used when testing k groups for homogeneity in response
rates. The null hypothesis is r₁=...=r_(k), and the alternative is that
there exists a pair of groups with differing rates. `ZSquared`
implements the square of a normally distributed random variable with
mean \\\mu\\ and standard deviation \\\sigma^2\\.

## Usage

``` r
ChiSquared(df)

# S4 method for class 'ChiSquared'
quantile(x, probs, n, theta, ...)

# S4 method for class 'ChiSquared,numeric'
simulate(object, nsim, n, theta, seed = NULL, ...)
```

## Arguments

- df:

  number of degrees of freedom

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

  object of class `ChiSquared`

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
datadist <- ChiSquared(df=4)
```
