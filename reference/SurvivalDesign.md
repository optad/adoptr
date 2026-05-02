# SurvivalDesign

`SurvivalDesign` is a function that converts an arbitrary design to a
survival design.

## Usage

``` r
SurvivalDesign(design, event_rate)

# S4 method for class 'TwoStageDesign'
SurvivalDesign(design, event_rate)

# S4 method for class 'TwoStageDesign'
TwoStageDesign(n1, event_rate)

# S4 method for class 'OneStageDesign'
OneStageDesign(n, event_rate)

# S4 method for class 'OneStageDesign'
SurvivalDesign(design, event_rate)

# S4 method for class 'GroupSequentialDesign'
GroupSequentialDesign(n1, event_rate)

# S4 method for class 'GroupSequentialDesign'
SurvivalDesign(design, event_rate)
```

## Arguments

- design:

  design that should be converted to a survival design

- event_rate:

  probability that a subject in either group will eventually have an
  event

- n1:

  design object to convert (overloaded from `TwoStageDesign`)

- n:

  design object to convert (overloaded from `TwoStageDesign`)

## Value

Converts any type of design to a survival design

## Examples

``` r
design <- get_initial_design(0.4, 0.025, 0.1)
SurvivalDesign(design, 0.8)
#> TwoStageDesignSurvival<n_events1=72;0.0<=x1<=2.2;n_events2=33-182> 

design_os <- get_initial_design(0.4, 0.025, 0.1, type_design = "one-stage")
design_gs <- get_initial_design(0.4, 0.025, 0.1, type_design = "group-sequential")

OneStageDesign(design_os, 0.7)
#> OneStageDesignSurvival<n_events=131;c=1.96> 

GroupSequentialDesign(design_gs, 0.8)
#> GroupSequentialDesignSurvival<n_events1=72;0.0<=x1<=2.2;n_events2=72> 
```
