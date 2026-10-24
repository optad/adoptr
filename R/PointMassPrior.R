#' Univariate discrete point mass priors
#'
#' \code{PointMassPrior} is a sub-class of \code{\link[=Prior-class]{Prior}}
#' representing a univariate prior over a discrete set of points with positive
#' probability mass.
#'
#' @slot theta cf. parameter 'theta'
#' @slot mass cf. parameter 'mass'
#'
#' @seealso To represent continuous prior distributions use \code{\link{ContinuousPrior}}.
#'
#' @aliases PointMassPrior
#' @exportClass PointMassPrior
setClass("PointMassPrior", representation(
        theta = "numeric",
        mass  = "numeric"
    ),
    contains = "Prior")



#' @param theta numeric vector of pivot points with positive prior mass
#' @param mass numeric vector of probability masses at the pivot points
#'     (must sum to 1)
#' @template label
#'
#' @return an object of class \code{PointMassPrior}, \code{theta} is
#' automatically sorted in ascending order
#'
#' @examples
#' PointMassPrior(c(0, .5), c(.3, .7))
#'
#' @rdname PointMassPrior-class
#' @export
PointMassPrior <- function(theta, mass, label = NA_character_) {
    if (length(theta) != length(mass))
        stop("theta and mass must be of same length")
    if (any(is.na(mass)) || any(mass < 0))
        stop("mass must be non-negative")
    if (round(sum(mass), 10) != 1)
        stop("mass must sum to one")
    ord <- order(theta)
    new("PointMassPrior", theta = theta[ord], mass = mass[ord],
        label = label)
}




setMethod("print", signature('PointMassPrior'), function(x, ...) {
    name <- if (!is.na(x@label)) x@label else 'PointMass'
    if (length(x@theta) == 1)
        glue::glue("{name}<{sprintf('%.2f',x@theta)}>")
    else
        paste0(
            glue::glue("{name}<"),
            paste0(glue::glue("Pr[{sprintf('%.2f',x@theta)}]={sprintf('%3.2f',x@mass)}"), collapse = ";"),
            ">",
            collapse = ""
        )
})





#' @examples
#' bounds(PointMassPrior(c(0, .5), c(.3, .7)))
#' # > 0.0 0.5
#'
#' @rdname bounds
#' @export
setMethod("bounds", signature("PointMassPrior"),
    function(dist, ...) range(dist@theta))


#' @examples
#' expectation(PointMassPrior(c(0, .5), c(.3, .7)), identity)
#' # > .35
#'
#' @rdname expectation
#' @export
setMethod("expectation", signature("PointMassPrior", "function"),
    function(dist, f, ...) sum(dist@mass * sapply(dist@theta, f, ...)) )


#' @examples
#' tmp <- condition(PointMassPrior(c(0, .5), c(.3, .7)), c(-1, .25))
#' expectation(tmp, identity) # 0
#'
#' @rdname condition
#' @export
setMethod("condition", signature("PointMassPrior", "numeric"),
    function(dist, interval, ...) {
        if (length(interval) != 2)
            stop("interval must be of length 2")
        if (any(!is.finite(interval)))
            stop("interval must be finite")
        if (diff(interval) < 0)
            stop("interval[2] must be larger or equal to interval[1]")
        epsilon <- sqrt(.Machine$double.eps)
        # find indices of pivots wihtin interval (up to machine precision!)
        idx <- (interval[1] - dist@theta <= epsilon) & (dist@theta - interval[2] <= epsilon)
        if (!any(idx) || sum(dist@mass[idx]) <= 0)
            stop("the prior has no mass in the given interval")
        # re-normalize and return
        return(PointMassPrior(
            dist@theta[idx], dist@mass[idx] / sum(dist@mass[idx])
        ))
    })


#' @examples
#' predictive_pdf(Normal(), PointMassPrior(.3, 1), 1.5, 20) # ~.343
#'
#' @rdname predictive_pdf
#' @export
setMethod("predictive_pdf", signature("DataDistribution", "PointMassPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        grid <- grid_x_theta(x1, prior@theta)
        (matrix(probability_density_function(dist, grid$x, n1, grid$theta), nrow = length(x1)) %*% prior@mass)[, 1]
    })


#' @examples
#' predictive_cdf(Normal(), PointMassPrior(.0, 1), 0, 20) # .5
#'
#' @rdname predictive_cdf
#' @export
setMethod("predictive_cdf", signature("DataDistribution", "PointMassPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        grid <- grid_x_theta(x1, prior@theta)
        (matrix(cumulative_distribution_function(dist, grid$x, n1, grid$theta), nrow = length(x1)) %*% prior@mass)[, 1]
    })


#' @examples
#' posterior(Normal(), PointMassPrior(0, 1), 2, 20)
#'
#' @rdname posterior
#' @export
setMethod("posterior", signature("DataDistribution", "PointMassPrior", "numeric"),
    function(dist, prior, x1, n1, ...) {
        if (length(prior@theta) == 1) return(prior) # shortcut
        mass  <- prior@mass * probability_density_function(dist, x1, n1, prior@theta)
        total <- sum(mass)
        if (!is.finite(total) || total <= 0)
            stop("posterior is undefined: the observation has zero likelihood under every support point of the prior")
        PointMassPrior(prior@theta, mass / total) # normalize
    })
