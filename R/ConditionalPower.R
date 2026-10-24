#' (Conditional) Power of a Design
#'
#' This score evaluates \ifelse{html}{\out{P[X<sub>2</sub> > c2(design, X<sub>1</sub>) | X<sub>1</sub> = x<sub>1</sub>]}}{\eqn{\boldsymbol{P}[X_2 > c_2(design, X_1)|X_1 = x_1]}}.
#' Note that the distribution of \ifelse{html}{\out{X<sub>2</sub>}}{\eqn{X_2}} is the posterior predictive after
#' observing \ifelse{html}{\out{X<sub>1</sub> = x<sub>1</sub>}}{\eqn{X_1 = x_1}}.
#'
#' @template dist
#' @template prior
#' @template design
#' @template s
#' @template x1
#' @template optimization
#' @template dotdotdot
#'
#' @seealso \code{\link{Scores}}
#'
#' @examples
#' prior <- PointMassPrior(.4, 1)
#' cp <- ConditionalPower(Normal(), prior)
#' evaluate(
#'    cp,
#'    TwoStageDesign(50, .0, 2.0, 50, 2.0, order = 5L),
#'    x1 = 1
#' )
#' # these two are equivalent:
#' expected(cp, Normal(), prior)
#' Power(Normal(), prior)
#'
#' @aliases ConditionalPower
#' @exportClass ConditionalPower
setClass("ConditionalPower", representation(
        distribution = "DataDistribution",
        prior        = "Prior"
    ),
    contains = "ConditionalScore")


#' @template label
#'
#' @rdname ConditionalPower-class
#' @export
ConditionalPower <- function(dist, prior, label = "Pr[x2>=c2(x1)|x1]") {
    new("ConditionalPower", distribution = dist, prior = prior, label = label)
}



#' @rdname ConditionalPower-class
#' @export
Power <- function(dist, prior, label = "Pr[x2>=c2(x1)]") {
    expected(ConditionalPower(dist, prior), dist, prior, label = label)
}



# internal: conditional power for a discrete prior, vectorised over x1.
#
# For every stage-one outcome x1[i] the posterior over the support points of the
# prior is computed and the conditional power
#   sum_j posterior_ij * (1 - F(c2[i] | n2[i], theta_j))
# is returned. The two data-distribution calls (pdf for the posterior, cdf for
# the power) are made once on the whole (x1 x theta) grid instead of once per
# x1 (and per support point), which is considerably faster inside the
# optimization loop. The data distribution methods must be vectorised in x, n
# (n1 is scalar, n2 may be a vector of the same length as x1) and theta.
conditional_power_point_mass <- function(dist, prior, x1, n1, c2, n2) {
    m     <- length(x1)
    k     <- length(prior@theta)
    idx_x <- rep(seq_len(m), times = k)
    idx_t <- rep(seq_len(k), each = m)
    theta <- prior@theta[idx_t]

    if (k == 1L) {
        post <- matrix(1, nrow = m, ncol = 1L)
    } else {
        dens  <- matrix(probability_density_function(dist, x1[idx_x], n1, theta), nrow = m)
        post  <- dens * rep(prior@mass, each = m)
        total <- rowSums(post)
        if (any(!is.finite(total) | total <= 0))
            stop("posterior is undefined: an observation has zero likelihood under every support point of the prior")
        post  <- post / total
    }

    n2_grid <- if (length(n2) == 1L) n2 else n2[idx_x]
    exceed  <- 1 - matrix(
        cumulative_distribution_function(dist, c2[idx_x], n2_grid, theta),
        nrow = m
    )
    return(rowSums(post * exceed))
}



#' @rdname ConditionalPower-class
#' @export
setMethod("evaluate", signature("ConditionalPower", "TwoStageDesign"),
    function(s, design, x1, optimization = FALSE, ...) {
        point_mass <- is(s@prior, "PointMassPrior")
        if (optimization) { # use design-specific implementation
            pivots <- scaled_integration_pivots(design)
            idx    <- vapply(x1, function(x) which.min(abs(x - pivots)), integer(1)) # more robust than ==; know x1 are pivots!
            n2_piv <- if (is(design, 'GroupSequentialDesign')) design@n2_pivots else design@n2_pivots[idx]
            if (point_mass) {
                res <- numeric(length(x1)) # x1 < c1f: 0
                res[x1 > design@c1e] <- 1
                cont <- which(x1 >= design@c1f & x1 <= design@c1e)
                if (length(cont) > 0) {
                    res[cont] <- conditional_power_point_mass(
                        s@distribution, s@prior,
                        pivots[idx[cont]], design@n1,
                        design@c2_pivots[idx[cont]],
                        if (length(n2_piv) == 1L) n2_piv else n2_piv[cont]
                    )
                }
                return(res)
            }
            sapply(
                seq_along(x1),
                function(i) {
                    if (x1[i] < design@c1f) return(0)
                    if (x1[i] > design@c1e) return(1)
                    return(
                        expectation(
                            posterior( # candidate for memoisation
                                s@distribution,
                                s@prior,
                                pivots[idx[i]],
                                design@n1,
                                ...
                            ),
                            function(theta) 1 - cumulative_distribution_function( # candidate for memoisation
                                s@distribution,
                                design@c2_pivots[idx[i]],
                                if (is(design, 'GroupSequentialDesign')) design@n2_pivots else design@n2_pivots[idx[i]],
                                theta
                            )
                        )
                    )
                }
            )
        } else {
            if (point_mass) {
                return(conditional_power_point_mass(
                    s@distribution, s@prior,
                    x1, n1(design, round = TRUE),
                    c2(design, x1), n2(design, x1, round = TRUE)
                ))
            }
            sapply(
                x1,
                function(x1) expectation(
                    posterior(s@distribution, s@prior, x1, n1(design, round = TRUE), ...),
                    function(theta) 1 - cumulative_distribution_function(s@distribution, c2(design, x1), n2(design, x1, round = TRUE), theta)
                )
            )
        }
    })



# not user facing!
setMethod(".evaluate", signature("ConditionalPower", "TwoStageDesign"),
    function(s, design, ...) {
        pivots <- scaled_integration_pivots(design)
        n2_piv <- design@n2_pivots # length one for group-sequential designs
        if (is(s@prior, "PointMassPrior")) {
            continuation <- conditional_power_point_mass(
                s@distribution, s@prior,
                pivots, design@n1, design@c2_pivots, n2_piv
            )
        } else {
            n2 <- if (is(design, 'GroupSequentialDesign')) function(i) design@n2_pivots else function(i) design@n2_pivots[i]
            continuation <- sapply(
                seq_along(pivots),
                function(i) {
                    expectation(
                        posterior(s@distribution, s@prior, pivots[i], design@n1),
                        function(theta) 1 - cumulative_distribution_function(s@distribution, design@c2_pivots[i], n2(i), theta)
                    )
                }
            )
        }
        list(
            early_futility = 0,
            pivots         = continuation,
            early_efficacy = 1
        )
    })

setMethod(".evaluate", signature("ConditionalPower", "OneStageDesign"),
    function(s, design, ...) {
        list(
            early_futility = 0,
            pivots         = NaN,
            early_efficacy = 1
        )
    })
