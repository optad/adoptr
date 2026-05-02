# F-Distribution

Implements the F-distribution used for an ANOVA or for the comparison of
the fit of two nested regression models. In both cases, the test
statistic follows a F-distribution. `NestedModel` is used to compare the
fit of two regression models, where one model contains the independent
variables of the smaller model as a subset. Then, one can use ANOVA to
determine whether more variance can be explained by adding more
independent variables. In the class `ANOVA`, the number of independent
variables of the smaller model is set to \\1\\ in order to match the
degrees of freedom and we obtain a one-way ANOVA.

## Usage

``` r
NestedModels(p_inner, p_outer)

# S4 method for class 'NestedModels'
quantile(x, probs, n, theta, ...)

# S4 method for class 'NestedModels,numeric'
simulate(object, nsim, n, theta, seed = NULL, ...)
```

## Arguments

- p_inner:

  number of independent variables in smaller model

- p_outer:

  number of independent variables in bigger model

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

  object of class `NestedModels`

- nsim:

  number of simulation runs

- seed:

  random seed

## Slots

- `p_inner`:

  number of parameters in smaller model

- `p_outer`:

  number of parameters in bigger model

## See also

See
[`probability_density_function`](https://optad.github.io/adoptr/reference/probability_density_function.md)
and
[`cumulative_distribution_function`](https://optad.github.io/adoptr/reference/cumulative_distribution_function.md)
to evaluate the pdf and the cdf, respectively. Use
[`ANOVA`](https://optad.github.io/adoptr/reference/ANOVA-class.md) for
detailed information of ANOVA.

## Examples

``` r
model <- NestedModels(2, 4)
```
