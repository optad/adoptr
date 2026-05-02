# Two-stage design for time-to-event-endpoints

When conducting a study with time-to-event endpoints, the main interest
is not the sample size, but the number of overall necessary events.
Thus, [adoptr](https://optad.github.io/adoptr/reference/adoptr.md) does
not use the sample size for calculating the design. Instead, it uses the
number of events directly. In the framework of
[adoptr](https://optad.github.io/adoptr/reference/adoptr.md), all the
calculations are done group-wise, where both of the groups are
equal-sized. This means, that the number of events
[adoptr](https://optad.github.io/adoptr/reference/adoptr.md) has
computed is only half of the overall number of necessary events. In
order to facilitate this issue, the look of the `summary` and `show`
functions have been changed in the survival analysis setting. The sample
size is implicitly determined by dividing the number of events by the
event rate. Survival objects are only created, when the argument
`event_rate` is not missing.

## Slots

- `event_rate`:

  probability that a subject in either group will eventually have an
  event

## See also

[`TwoStageDesign`](https://optad.github.io/adoptr/reference/TwoStageDesign-class.md)
for superclass and inherited methods
