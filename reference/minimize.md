# Find optimal two-stage design by constraint minimization

`minimize` takes an unconditional score and a constraint set (or no
constraint) and solves the corresponding minimization problem using
[`nloptr`](https://cran.r-project.org/package=nloptr) (using COBYLA by
default). An initial design has to be defined. It is also possible to
define lower- and upper-boundary designs. If this is not done, the
boundaries are determined automatically heuristically.

## Usage

``` r
minimize(
  objective,
  subject_to,
  initial_design,
  lower_boundary_design = get_lower_boundary_design(initial_design),
  upper_boundary_design = get_upper_boundary_design(initial_design),
  c2_decreasing = FALSE,
  check_constraints = TRUE,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
  ...
)
```

## Arguments

- objective:

  objective function

- subject_to:

  constraint collection

- initial_design:

  initial guess (x0 for nloptr)

- lower_boundary_design:

  design specifying the lower boundary.

- upper_boundary_design:

  design specifying the upper boundary

- c2_decreasing:

  if TRUE, the c2_pivots are forced to be monotonically decreasing

- check_constraints:

  if TRUE, it is checked if constrains are fulfilled

- opts:

  options list passed to nloptr

- ...:

  further optional arguments passed to
  [`nloptr`](https://astamm.github.io/nloptr/reference/nloptr.html)

## Value

a list with elements:

- design:

  The resulting optimal design

- nloptr_return:

  Output of the corresponding nloptr call

- call_args:

  The arguments given to the optimization call

## Examples

``` r
# Define Type one error rate
toer <- Power(Normal(), PointMassPrior(0.0, 1))

# Define Power at delta = 0.4
pow <- Power(Normal(), PointMassPrior(0.4, 1))

# Define expected sample size at delta = 0.4
ess <- ExpectedSampleSize(Normal(), PointMassPrior(0.4, 1))

# Compute design minimizing ess subject to power and toer constraints
# \donttest{
minimize(

   ess,

   subject_to(
      toer <= 0.025,
      pow  >= 0.9
   ),

   initial_design = TwoStageDesign(50, .0, 2.0, 60.0, 2.0, 5L)

)
#> $design
#> TwoStageDesign<n1=68;0.3<=x1<=2.3:n2=26-131> 
#> 
#> $nloptr_return
#> 
#> Call:
#> nloptr::nloptr(x0 = tunable_parameters(initial_design), eval_f = f_obj, 
#>     lb = tunable_parameters(lower_boundary_design), ub = tunable_parameters(upper_boundary_design), 
#>     eval_g_ineq = g_cnstr, opts = opts)
#> 
#> 
#> Minimization using NLopt version 2.10.0 
#> 
#> NLopt solver status: 4 ( NLOPT_XTOL_REACHED: Optimization stopped because 
#> xtol_rel or xtol_abs (above) was reached. )
#> 
#> Number of Iterations....: 3990 
#> Termination conditions:  xtol_rel: 1e-05 maxeval: 10000 
#> Number of inequality constraints:  3 
#> Number of equality constraints:    0 
#> Optimal value of objective function:  99.2099508305843 
#> Optimal value of controls: 67.77852 0.2812201 2.26569 126.9927 111.6292 86.5138 57.48934 32.48944 2.668702 
#> 2.35649 1.822003 1.113608 0.3340881
#> 
#> 
#> 
#> $call_args
#> $call_args$objective
#> E[n(x1)] 
#> 
#> $call_args$subject_to
#> An object of class "ConstraintsCollection"
#> Slot "unconditional_constraints":
#> [[1]]
#> Pr[x2>=c2(x1)] <= 0.025 
#> 
#> [[2]]
#> -Pr[x2>=c2(x1)]  <= -0.9 
#> 
#> 
#> Slot "conditional_constraints":
#> list()
#> 
#> 
#> $call_args$initial_design
#> TwoStageDesign<n1=50;0.0<=x1<=2.0:n2=60> 
#> 
#> $call_args$lower_boundary_design
#> TwoStageDesign<n1=1;-2.0<=x1<=0.0:n2=1> 
#> 
#> $call_args$upper_boundary_design
#> TwoStageDesign<n1=250;2.0<=x1<=4.0:n2=300> 
#> 
#> $call_args$c2_decreasing
#> [1] FALSE
#> 
#> $call_args$check_constraints
#> [1] TRUE
#> 
#> $call_args$opts
#> $call_args$opts$algorithm
#> [1] "NLOPT_LN_COBYLA"
#> 
#> $call_args$opts$xtol_rel
#> [1] 1e-05
#> 
#> $call_args$opts$maxeval
#> [1] 10000
#> 
#> 
#> 
#> attr(,"class")
#> [1] "adoptrOptimizationResult" "list"                    
# }
```
