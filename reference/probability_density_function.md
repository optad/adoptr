# Probability density function

`probability_density_function` evaluates the probability density
function of a specific distribution `dist` at a point `x`.

## Usage

``` r
probability_density_function(dist, x, n, theta, ...)

# S4 method for class 'Binomial,numeric,numeric,numeric'
probability_density_function(dist, x, n, theta, ...)

# S4 method for class 'ChiSquared,numeric,numeric,numeric'
probability_density_function(dist, x, n, theta, ...)

# S4 method for class 'NestedModels,numeric,numeric,numeric'
probability_density_function(dist, x, n, theta, ...)

# S4 method for class 'Normal,numeric,numeric,numeric'
probability_density_function(dist, x, n, theta, ...)

# S4 method for class 'Student,numeric,numeric,numeric'
probability_density_function(dist, x, n, theta, ...)

# S4 method for class 'Survival,numeric,numeric,numeric'
probability_density_function(dist, x, n, theta, ...)
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

value of the probability density function at point `x`.

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
probability_density_function(Binomial(.2, FALSE), 1, 50, .3)
#> [1] 0.0008519612

probability_density_function(Pearson2xK(3), 1, 30, get_tau_Pearson2xK(c(0.3, 0.4, 0.7, 0.2)))
#> [1] 0.003505548
probability_density_function(ZSquared(TRUE), 1, 35, get_tau_ZSquared(0.4, 1))
#> [1] 0.1646112


probability_density_function(ANOVA(3), 1, 30, get_tau_ANOVA(c(0.3, 0.4, 0.7, 0.2)))
#> [1] 0.2513264

probability_density_function(Normal(), 1, 50, .3)
#> [1] 0.3520653

probability_density_function(Student(TRUE), 1, 40, 1.1)
#> [1] 0.0001946335

probability_density_function(Survival(0.6,TRUE),0.75,50,0.9)
#> [1] 0.1765677
```
