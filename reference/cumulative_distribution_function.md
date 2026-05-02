# Cumulative distribution function

`cumulative_distribution_function` evaluates the cumulative distribution
function of a specific distribution `dist` at a point `x`.

## Usage

``` r
cumulative_distribution_function(dist, x, n, theta, ...)

# S4 method for class 'Binomial,numeric,numeric,numeric'
cumulative_distribution_function(dist, x, n, theta, ...)

# S4 method for class 'ChiSquared,numeric,numeric,numeric'
cumulative_distribution_function(dist, x, n, theta, ...)

# S4 method for class 'NestedModels,numeric,numeric,numeric'
cumulative_distribution_function(dist, x, n, theta, ...)

# S4 method for class 'Normal,numeric,numeric,numeric'
cumulative_distribution_function(dist, x, n, theta, ...)

# S4 method for class 'Student,numeric,numeric,numeric'
cumulative_distribution_function(dist, x, n, theta, ...)

# S4 method for class 'Survival,numeric,numeric,numeric'
cumulative_distribution_function(dist, x, n, theta, ...)
```

## Arguments

- dist:

  a univariate
  [`distribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md)
  object

- x:

  outcome

- n:

  sample size

- theta:

  distribution parameter

- ...:

  further optional arguments

## Value

value of the cumulative distribution function at point `x`.

## Details

If the distribution is
[`Binomial`](https://optad.github.io/adoptr/reference/BinomialDataDistribution-class.md),
theta denotes the rate difference between intervention and control
group. Then, the mean is assumed to be √ n theta.

If the distribution is
[`Normal`](https://optad.github.io/adoptr/reference/NormalDataDistribution-class.md),
then the mean is assumed to be √ n theta.

## Examples

``` r
cumulative_distribution_function(Binomial(.1, TRUE), 1, 50, .3)
#> [1] 0.004310344

cumulative_distribution_function(Pearson2xK(3), 1, 30, get_tau_Pearson2xK(c(0.3,0.4,0.7,0.2)))
#> [1] 0.001966853
cumulative_distribution_function(ZSquared(TRUE), 1, 35, get_tau_ZSquared(0.4, 1))
#> [1] 0.2466166


cumulative_distribution_function(ANOVA(3), 1, 30, get_tau_ANOVA(c(0.3, 0.4, 0.7, 0.2)))
#> [1] 0.2402678

cumulative_distribution_function(Normal(), 1, 50, .3)
#> [1] 0.3085375

cumulative_distribution_function(Student(two_armed = FALSE), .75, 50, .9)
#> [1] 1.062003e-08

cumulative_distribution_function(Survival(0.6,TRUE),0.75,50,0.9)
#> [1] 0.899164
```
