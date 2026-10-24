#' F-Distribution
#'
#' Implements the F-distribution used for an ANOVA or for the comparison of the fit of two
#' nested regression models. In both cases, the test statistic follows a F-distribution.
#' \code{NestedModel} is used to compare the fit of two regression models, where one model contains
#' the independent variables of the smaller model as a subset. Then, one can use ANOVA to determine
#' whether more variance can be explained by adding more independent variables.
#' In the class \code{ANOVA}, the number of independent variables of the smaller model is set to \eqn{1}
#' in order to match the degrees of freedom and we obtain a one-way ANOVA.
#'
#' @slot p_inner number of parameters in smaller model
#' @slot p_outer number of parameters in bigger model
#'
#' @template DataDistributionTemplate
#'
#' @include DataDistribution.R
#'
#' @rdname NestedModels-class
#' @exportClass NestedModels
setClass("NestedModels", representation(
  p_inner  = "numeric",
  p_outer  =  "numeric"),
  contains = "DataDistribution")


#' @include DataDistribution.R
#'
#' @rdname ANOVA-class
#' @exportClass ANOVA
setClass("ANOVA", contains = "NestedModels")



#' @param p_inner number of independent variables in smaller model
#' @param p_outer number of independent variables in bigger model
#'
#' @examples
#' model <- NestedModels(2, 4)
#'
#' @seealso See \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively. Use \code{\link{ANOVA}} for detailed information
#'    of ANOVA.
#'
#' @rdname NestedModels-class
#' @export
NestedModels <- function(p_inner, p_outer) {
  tol <- sqrt(.Machine$double.eps)
  if (any(p_inner < 0 | abs(p_inner - round(p_inner)) > tol,
          p_outer < 0 | abs(p_outer - round(p_outer)) > tol)) {
    stop("The numbers of groups must be natural numbers.")
  }
  if (any(p_outer <= p_inner))
    stop("p_outer must be larger than p_inner.")

  new("NestedModels", p_inner = p_inner, p_outer = p_outer)
}


#' Analysis of Variance
#'
#' ANOVA is used to test whether there is a significant difference between the means of groups.
#' The sample size which \code{adoptr} returns is the group wise sample size.
#' The function \code{get_tau_ANOVA} is used to obtain a parameter \eqn{\tau},
#' which is used in the same way as \eqn{\theta} to describe the difference of
#' means between the groups.
#'
#' @param n_groups number of groups to be compared
#'
#' @examples
#' model <- ANOVA(3L)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively. Use \code{\link{NestedModels}} to get insights
#'    in the implementation of \code{ANOVA}.
#'
#' @rdname ANOVA-class
#' @export
ANOVA <- function(n_groups) {
  if (n_groups < 2 || abs(n_groups - round(n_groups)) > sqrt(.Machine$double.eps))
    stop("The number of groups must be a natural number of at least two.")
  new("ANOVA", p_outer = n_groups, p_inner = 1L)
}

#' @param means vector denoting the mean per group
#' @param common_sd standard deviation of the groups
#'
#' @examples
#' H1 <- PointMassPrior(get_tau_ANOVA(c(0.4, 0.8, 0.5)), 1)
#'
#' @rdname ANOVA-class
#'
#' @export
get_tau_ANOVA <- function(means, common_sd = 1) {
  mean((means - mean(means))^2) / common_sd^2
}


# internal: arguments of the (non-central) F-distribution, recycled to a
# common length. The sample size is converted to the overall sample size.
# Vectorised in x, n and theta (e.g. needed when a posterior over several
# support points is evaluated at a single x, or when n varies across
# simulation runs).
nested_models_args <- function(dist, x, n, theta) {
  len <- max(length(x), length(n), length(theta))
  n   <- rep_len(n, len) * (dist@p_outer - dist@p_inner + 1)
  list(
    x     = rep_len(x, len),
    n     = n,
    theta = rep_len(theta, len),
    df1   = dist@p_outer - dist@p_inner,
    df2   = pmax(1, n - dist@p_outer)
  )
}


#' @examples
#' probability_density_function(ANOVA(3), 1, 30, get_tau_ANOVA(c(0.3, 0.4, 0.7, 0.2)))
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function", signature("NestedModels", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
            a   <- nested_models_args(dist, x, n, theta)
            ret <- numeric(length(a$x)) # density is zero at +-Inf
            fin <- !is.infinite(a$x)
            ret[fin] <- stats::df(a$x[fin],
                                  df1 = a$df1,
                                  df2 = a$df2[fin],
                                  ncp = a$n[fin] * a$theta[fin])
            return(ret)
          })


#' @examples
#' cumulative_distribution_function(ANOVA(3), 1, 30, get_tau_ANOVA(c(0.3, 0.4, 0.7, 0.2)))
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function", signature("NestedModels", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
            a   <- nested_models_args(dist, x, n, theta)
            ret <- numeric(length(a$x))
            ret[which(a$x == Inf)] <- 1
            fin <- !is.infinite(a$x)
            ret[fin] <- stats::pf(a$x[fin],
                                  df1 = a$df1,
                                  df2 = a$df2[fin],
                                  ncp = a$n[fin] * a$theta[fin])
            return(ret)
          })


#' @param probs vector of probabilities
#'
#' @rdname NestedModels-class
#' @export
setMethod("quantile", signature("NestedModels"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
            n <- n * (x@p_outer - x@p_inner + 1)
            return(stats::qf(probs, df1 = x@p_outer - x@p_inner,
                             df2 = pmax(1, n - x@p_outer), ncp = n * theta))
          })



#' @rdname NestedModels-class
#'
#' @param object object of class \code{NestedModels}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("NestedModels", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
            if (!is.null(seed)) set.seed(seed)
            n <- n * (object@p_outer - object@p_inner + 1)
            return(stats::rf(nsim, df1 = object@p_outer - object@p_inner,
                             df2 = pmax(1, n - object@p_outer),  ncp = n * theta))
          })

setMethod("print", signature('NestedModels'), function(x, ...) {
  glue::glue(
    "{class(x)[1]}<p_inner={x@p_inner},p_outer={x@p_outer}>"
  )
})

setMethod("print", signature('ANOVA'), function(x, ...) {
  glue::glue(
    "{class(x)[1]}<n_groups={x@p_outer}>"
  )
})
