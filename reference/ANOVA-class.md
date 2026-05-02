# Analysis of Variance

ANOVA is used to test whether there is a significant difference between
the means of groups. The sample size which `adoptr` returns is the group
wise sample size. The function `get_tau_ANOVA` is used to obtain a
parameter \\\tau\\, which is used in the same way as \\\theta\\ to
describe the difference of means between the groups.

## Usage

``` r
ANOVA(n_groups)

get_tau_ANOVA(means, common_sd = 1)
```

## Arguments

- n_groups:

  number of groups to be compared

- means:

  vector denoting the mean per group

- common_sd:

  standard deviation of the groups

## See also

see
[`probability_density_function`](https://optad.github.io/adoptr/reference/probability_density_function.md)
and
[`cumulative_distribution_function`](https://optad.github.io/adoptr/reference/cumulative_distribution_function.md)
to evaluate the pdf and the cdf, respectively. Use
[`NestedModels`](https://optad.github.io/adoptr/reference/NestedModels-class.md)
to get insights in the implementation of `ANOVA`.

## Examples

``` r
model <- ANOVA(3L)

H1 <- PointMassPrior(get_tau_ANOVA(c(0.4, 0.8, 0.5)), 1)
```
