# Adaptive Optimal Two-Stage Designs

The adoptr package provides functionality to explore custom optimal
two-stage designs for one- or two-arm superiority tests. For more
details on the theoretical background see <doi:10.1002/sim.8291> and
<doi:10.18637/jss.v098.i09>. adoptr makes heavy use of the S4 class
system. A good place to start learning about it can be found
[here](http://adv-r.had.co.nz/OO-essentials.md).

## Quickstart

For a sample workflow and a quick demo of the capabilities, see
[here](https://optad.github.io/adoptr/articles/adoptr.html).

A more detailed description of the background and the usage of adoptr
can be found
[here](https://optad.github.io/adoptr/articles/adoptr_jss.html) or here
<doi:10.18637/jss.v098.i09> .

A variety of examples is presented in the validation report hosted
[here](https://optad.github.io/adoptr-validation-report/).

## Designs

adoptr currently supports
[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md),
[`GroupSequentialDesign`](https://optad.github.io/adoptr/reference/GroupSequentialDesign-class.md),
and
[`OneStageDesign`](https://optad.github.io/adoptr/reference/OneStageDesign-class.md).

## Data distributions

The implemented data distributions are
[`Normal`](https://optad.github.io/adoptr/reference/NormalDataDistribution-class.md),
[`Binomial`](https://optad.github.io/adoptr/reference/BinomialDataDistribution-class.md),
[`Student`](https://optad.github.io/adoptr/reference/StudentDataDistribution-class.md),
[`Survival`](https://optad.github.io/adoptr/reference/SurvivalDataDistribution-class.md),
[`ChiSquared`](https://optad.github.io/adoptr/reference/ChiSquaredDataDistribution-class.md)
(including
[`Pearson2xK`](https://optad.github.io/adoptr/reference/Pearson2xK-class.md)
and
[`ZSquared`](https://optad.github.io/adoptr/reference/ZSquared-class.md))
and [`ANOVA`](https://optad.github.io/adoptr/reference/ANOVA-class.md).

## Priors

Both
[`ContinuousPrior`](https://optad.github.io/adoptr/reference/ContinuousPrior-class.md)
and
[`PointMassPrior`](https://optad.github.io/adoptr/reference/PointMassPrior-class.md)
are supported for the single parameter of a
[`DataDistribution`](https://optad.github.io/adoptr/reference/DataDistribution-class.md).

## Scores

See [`Scores`](https://optad.github.io/adoptr/reference/Scores.md) for
information on the basic system of representing scores. Available scores
are
[`ConditionalPower`](https://optad.github.io/adoptr/reference/ConditionalPower-class.md),
[`ConditionalSampleSize`](https://optad.github.io/adoptr/reference/ConditionalSampleSize-class.md),
[`Power`](https://optad.github.io/adoptr/reference/ConditionalPower-class.md),
and
[`ExpectedSampleSize`](https://optad.github.io/adoptr/reference/ConditionalSampleSize-class.md).

## See also

Useful links:

- <https://github.com/optad/adoptr>

- <https://optad.github.io/adoptr/>

- Report bugs at <https://github.com/optad/adoptr/issues>

## Author

**Maintainer**: Maximilian Pilz <maximilian.pilz@th-nuernberg.de>
([ORCID](https://orcid.org/0000-0002-9685-1613))

Authors:

- Maximilian Pilz <maximilian.pilz@th-nuernberg.de>
  ([ORCID](https://orcid.org/0000-0002-9685-1613))

- Kevin Kunzmann <kevin.kunzmann@boehringer-ingelheim.com>
  ([ORCID](https://orcid.org/0000-0002-1140-7143)) \[copyright holder\]

- Jan Meis <meis@imbi.uni-heidelberg.de>
  ([ORCID](https://orcid.org/0000-0001-5407-7220))

- Nico Bruder <bruder@imbi.uni-heidelberg.de>
  ([ORCID](https://orcid.org/0009-0004-9522-2075))
