#' Chi-Squared data distribution
#'
#' Implements a chi-squared distribution. The classes \code{Pearson2xk}
#' and \code{ZSquared} are subclasses, used in two different situations.
#' \code{Pearson2xK} is used when testing k groups for homogeneity in
#' response rates. The null hypothesis is
#' \ifelse{html}{\out{r<sub>1</sub>}=...=\out{r<sub>k</sub>}}{\eqn{r_1=...=r_k}}, and the
#' alternative is that there exists a pair of groups with differing rates.
#' \code{ZSquared} implements the square of a normally distributed random variable
#' with mean \eqn{\mu} and standard deviation \eqn{\sigma^2}.
#'
#' @template DataDistributionTemplate
#'
#' @rdname ChiSquaredDataDistribution-class
#' @exportClass ChiSquared
setClass("ChiSquared", representation(
  df  = "numeric",
  multiplier = "numeric"),
  contains = "DataDistribution")


#' Pearson's chi-squared test for contingency tables
#'
#' @include DataDistribution.R
#'
#' @rdname Pearson2xK-class
#' @exportClass Pearson2xK
setClass("Pearson2xK", contains = "ChiSquared")


#' Distribution class of a squared normal distribution
#'
#' @include DataDistribution.R
#'
#' @rdname ZSquared-class
#' @exportClass ZSquared
setClass("ZSquared", contains = "ChiSquared")


#' @param df number of degrees of freedom
#'
#' @examples
#' datadist <- ChiSquared(df=4)
#'
#' @seealso see \code{\link{probability_density_function}} and
#'    \code{\link{cumulative_distribution_function}} to evaluate the pdf
#'    and the cdf, respectively.
#'
#' @rdname ChiSquaredDataDistribution-class
#' @export
ChiSquared <- function(df) {
  if (df < 0 || abs(df - round(df)) > sqrt(.Machine$double.eps))
    stop("The degrees of freedom must be natural numbers.")
  new("ChiSquared", df = df, multiplier = 1)
}


#' Pearson's chi-squared test for 2 x k contingency tables
#'
#' When we test for homogeneity of rates in a k-armed trial with binary endpoints, 
#' the test statistic is chi-squared distributed with \eqn{k-1} degrees of 
#' freedom under the null. Under the alternative, the statistic is chi-squared 
#' distributed with a non-centrality parameter \eqn{\lambda}. 
#' The function \code{get_tau_Pearson2xk} then computes \eqn{\tau}, such that 
#' \eqn{\lambda} is given as \eqn{n \cdot \tau}, where \eqn{n} is the number of 
#' subjects per group. In \code{adoptr}, \eqn{\tau} is used in the same way as \eqn{\theta}
#' in the case of the normally distributed test statistic.
#'
#' @param n_groups number of groups considered for testing procedure
#'
#' @examples
#' pearson <- Pearson2xK(3)
#'
#'
#' @rdname Pearson2xK-class
#' @export
Pearson2xK <- function(n_groups) {
  if (n_groups < 0 || abs(n_groups - round(n_groups)) > sqrt(.Machine$double.eps))
    stop("The number of groups must be a natural number.")
  new("Pearson2xK", df = n_groups - 1L, multiplier = 1 / n_groups)
}

#' @param p_vector vector denoting the event rates per group
#'
#' @examples
#' H1 <- PointMassPrior(get_tau_Pearson2xK(c(.3, .25, .4)), 1)
#'
#' @rdname Pearson2xK-class
#'
#' @export
get_tau_Pearson2xK <- function(p_vector) {
  n_groups <- length(p_vector)
  mean_p <- mean(p_vector)
  deltas <- p_vector - mean_p
  tau <- (sum(deltas^2) / n_groups) / (mean_p * (1 - mean_p))
  tau
}


#' Squared normally distributed random variable
#'
#' Implementation of \eqn{Z^2}, where \eqn{Z} is normally distributed with mean 
#' \eqn{\mu} and variance \eqn{\sigma^2}. \eqn{Z^2} is chi-squared distributed 
#' with \eqn{1} degree of freedom and non-centrality parameter \eqn{(\mu/\sigma)^2}.
#' The function \code{get_tau_ZSquared} computes the factor \eqn{\tau=(\mu/\sigma)^2}, 
#' such that \eqn{\tau} is the equivalent of \eqn{\theta} in the normally 
#' distributed case. The square of a normal distribution \eqn{Z^2} can be used 
#' for two-sided hypothesis testing.
#'
#' @param two_armed logical indicating if a two-armed trial is regarded
#'
#' @examples
#' zsquared <- ZSquared(FALSE)
#'
#'
#' @rdname ZSquared-class
#' @export
ZSquared <- function(two_armed = TRUE) {
  new("ZSquared", df = 1L, multiplier = 1L + ifelse(two_armed, 1L, 0L))
}


#' @param mu mean of Z
#' @param sigma standard deviation of Z
#'
#' @examples
#' H1 <- PointMassPrior(get_tau_ZSquared(0.4, 1), 1)
#'
#' @rdname ZSquared-class
#'
#' @export
get_tau_ZSquared <- function(mu, sigma) {
  (mu / sigma)^2
}


#' @examples
#' probability_density_function(Pearson2xK(3), 1, 30, get_tau_Pearson2xK(c(0.3, 0.4, 0.7, 0.2)))
#' probability_density_function(ZSquared(4), 1, 35, get_tau_ZSquared(0.4))
#'
#'
#' @rdname probability_density_function
#' @export
setMethod("probability_density_function", signature("ChiSquared", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
            return(stats::dchisq(x, df = dist@df, ncp = n / dist@multiplier * theta))
          })


#' @examples
#' cumulative_distribution_function(Pearson2xK(3), 1, 30, get_tau_Pearson2xK(c(0.3,0.4,0.7,0.2)))
#' cumulative_distribution_function(ZSquared(4), 1, 35, get_tau_ZSquared(0.4))
#'
#'
#' @rdname cumulative_distribution_function
#' @export
setMethod("cumulative_distribution_function", signature("ChiSquared", "numeric", "numeric", "numeric"),
          function(dist, x, n, theta, ...) {
            return(stats::pchisq(x, df = dist@df, ncp = n / dist@multiplier * theta))
          })


#' @param probs vector of probabilities
#' @rdname ChiSquaredDataDistribution-class
#' @export
setMethod("quantile", signature("ChiSquared"),
          function(x, probs, n, theta, ...) { # must be x to conform with generic
            return(stats::qchisq(probs, df = x@df, ncp = n / x@multiplier * theta))
          })



#' @rdname ChiSquaredDataDistribution-class
#'
#' @param object object of class \code{ChiSquared}
#' @param nsim number of simulation runs
#' @param seed random seed
#'
#' @export
setMethod("simulate", signature("ChiSquared", "numeric"),
          function(object, nsim, n, theta, seed = NULL, ...) {
            if (!is.null(seed)) set.seed(seed)
            return(stats::rchisq(nsim, df = object@df, ncp = n / object@multiplier * theta))
          })

setMethod("print", signature('ChiSquared'), function(x, ...) {
  glue::glue(
    "{class(x)[1]}<df={x@df}>"
  )
})