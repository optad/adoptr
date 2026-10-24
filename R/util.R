# Compute weights an knots for Gauss-Legendre integration
#
# Implementation of Gauss Legendre integration rule and computation of integral
# based on implementation at https://rdrr.io/cran/SMR/src/R/GaussLegendre.R
# \code{GaussLegendreRule} gives the nodes and weights for the interval \eqn{(-1, 1)}.
# internal
GaussLegendreRule <- function(order) {
    order <- as.integer(order)
    if (order < 2) stop("At least two nodes are necessary for integration!")
    j   <- 1:(order - 1)
    mu0 <- 2
    b   <- j / (4 * j^2 - 1)^0.5
    A   <- rep(0, order * order)
    A[(order + 1) * (j - 1) + 2] <- b
    A[(order + 1) * j] <- b
    dim(A) <- c(order, order)
    sd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(sd$vectors[1, ]))
    w <- mu0 * w^2
    x <- rev(sd$values)
    return(data.frame(nodes = x, weights = w))
}



# integration via the Gauss-Legendre quadrature, internal
integrate_rule <- function(f, low, up, x, weights) {

  if (!(length(weights) == length(x)))
      stop("x and weights must be of same length")

  if (any(x < -1) | any(x > 1))
    stop("x must be in [-1, 1], is scaled automatically")

  if (any(weights <= 0))
    stop("weights must be positive")

  a  <- (up - low) / 2
  b  <- a + low
  ff <- f(a * x + b)

  return(a * sum(weights * ff))

}

# application of Gauss-Legendre quadrature on a set of pivots, internal
gauss_quad <- function(ff, low, up, weights) (up - low)/2 * sum(weights * ff)


# Fast replacement for expand.grid(x = x, theta = theta), internal
#
# Returns a list with the two columns of the grid (first factor varies
# fastest, exactly as in expand.grid). Avoids the construction of a data.frame,
# which is comparatively expensive and was called in every score evaluation
# during optimization.
grid_x_theta <- function(x, theta) {
    list(
        x     = rep(x, times = length(theta)),
        theta = rep(theta, each = length(x))
    )
}


setGeneric("as_character", function(x, ...) standardGeneric("as_character"))
