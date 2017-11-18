#' Accelerated univariate bisection method
#'
#' Find the root of a function between two points a and b using the
#' bisection method followed by the secant method.
#'
#' @param f function
#' @param a lower bound
#' @param b upper bound
#' @param tol tolerance, defaults to 10*.Machine$double.eps
#' @param maxit maximum number of iterations
#' @return a list
#' @export
#' @examples
#'
#' unibisect(function(x) x^2 - 2, 0, 2)
#' fast_unibisect(function(x) x^2 - 2, 0, 2)
#'
#'
#'
#' f <- function(x) pbeta(x, 90, 10) - .5
#' a <- 0; b <- 1
#' (out <- fast_unibisect(f, 0, 1))
#'
#' curve(f(x), col = "red", from = 0, to = 1)
#' with(out$evals, points(x, fx))
#'
#' library(ggplot2)
#' ggplot(out$evals, aes(x, fx, color = method)) +
#'   stat_function(fun = f, color = "black") +
#'   geom_point()
#'
#'
fast_unibisect <- function(f, a, b, tol = 10*.Machine$double.eps, bisections = 8L, maxit = 100L) {

  if(bisections >= maxit) return(unibisect(f, a, b, tol, maxit))

  # initialize x and fx
  x <- fx <- numeric(maxit)

  # check endpoints and return early if root there
  x[1] <- a;
  fa <- fx[1] <- f(a)
  if(abs(fa) <= tol) {
    return(list(
      root = a, val = fa,
      evals = data_frame(x = x[1], fx = fx[1], method = "bisection"),
      n_evals = 1
    ))
  }

  x[2] <- b
  fb <- fx[2] <- f(b)
  if(sign(fa) == sign(fb)) stop("f(a) and f(b) must have opposite signs.")
  if(abs(fb) <= tol) {
    return(list(
      root = b, val = fb,
      evals = data_frame(x = x[1:2], fx = fx[1:2], method = "bisection"),
      n_evals = 2
    ))
  }

  # compute first midpoint, return early if root found
  c  <- x[3] <- (a+b)/2
  fc <- fx[3] <- f(c)
  if(abs(fc) < tol) {
    return(list(
      root = c, val = fc,
      evals = data_frame(x = x[1:3], fx = fx[1:3], method = "bisection"),
      n_evals = 3
    ))
  }

  # bisect
  for(k in 4:bisections) {
    if (sign(fa) == sign(fc)) {
      a <- c; fa <- fc
    } else {
      b <- c; fb <- fc
    }
    c  <- x[k] <- (a+b)/2
    fc <- fx[k] <- f(c)
    n_evals <- k
    if(abs(fc) <= tol) break
  }

  # return if bisection successful
  if(abs(fc) < tol) {
    return(list(
      root = c, val = fc,
      evals = data_frame(x = x[1:n_evals], fx = fx[1:n_evals], method = "bisection"),
      n_evals = n_evals
    ))
  }

  # secant
  for(k in (bisections+1):maxit) {
    x[k]  <- x[k-1] - fx[k-1] / ((fx[k-1]-fx[k-2])/(x[k-1]-x[k-2]))
    fx[k] <- f(x[k])
    n_evals <- k
    if(abs(fx[k]) <= tol) break
  }


  # return
  list(
    root = x[k], f.root = fx[k],
    evals = data_frame(
      x = x[1:n_evals],
      fx = fx[1:n_evals],
      method = c(rep("bisection", bisections), rep("secant", n_evals-bisections))
    ),
    n_evals = n_evals
  )
}
