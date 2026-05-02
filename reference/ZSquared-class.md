# Distribution class of a squared normal distribution

Implementation of \\Z^2\\, where \\Z\\ is normally distributed with mean
\\\mu\\ and variance \\\sigma^2\\. \\Z^2\\ is chi-squared distributed
with \\1\\ degree of freedom and non-centrality parameter
\\(\mu/\sigma)^2\\. The function `get_tau_ZSquared` computes the factor
\\\tau=(\mu/\sigma)^2\\, such that \\\tau\\ is the equivalent of
\\\theta\\ in the normally distributed case. The square of a normal
distribution \\Z^2\\ can be used for two-sided hypothesis testing.

## Usage

``` r
ZSquared(two_armed = TRUE)

get_tau_ZSquared(mu, sigma)
```

## Arguments

- two_armed:

  logical indicating if a two-armed trial is regarded

- mu:

  mean of Z

- sigma:

  standard deviation of Z

## Examples

``` r
zsquared <- ZSquared(FALSE)


H1 <- PointMassPrior(get_tau_ZSquared(0.4, 1), 1)
```
