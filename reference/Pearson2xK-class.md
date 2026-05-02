# Pearson's chi-squared test for contingency tables

When we test for homogeneity of rates in a k-armed trial with binary
endpoints, the test statistic is chi-squared distributed with \\k-1\\
degrees of freedom under the null. Under the alternative, the statistic
is chi-squared distributed with a non-centrality parameter \\\lambda\\.
The function `get_tau_Pearson2xk` then computes \\\tau\\, such that
\\\lambda\\ is given as \\n \cdot \tau\\, where \\n\\ is the number of
subjects per group. In `adoptr`, \\\tau\\ is used in the same way as
\\\theta\\ in the case of the normally distributed test statistic.

## Usage

``` r
Pearson2xK(n_groups)

get_tau_Pearson2xK(p_vector)
```

## Arguments

- n_groups:

  number of groups considered for testing procedure

- p_vector:

  vector denoting the event rates per group

## Examples

``` r
pearson <- Pearson2xK(3)


H1 <- PointMassPrior(get_tau_Pearson2xK(c(.3, .25, .4)), 1)
```
