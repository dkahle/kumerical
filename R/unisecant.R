#' Univariate secant method
#'
#' Find the root of a function using the secant method.
#'
#' @param f function
#' @param x0 initial value
#' @param x1 second initial value
#' @param tol tolerance, defaults to 10*.Machine$double.eps
#' @param maxit maximum number of iterations
#' @return a list
#' @export
#' @examples
#'
#'
#' unisecant(function(x) x^2 - 2, 0, 1)
#'
#'
#' f  <- function(x) x^2 - 2
#' x0 <- 0; x1 <- 1
#' (out <- unisecant(f, x0, x1))
#'
#' curve(f(x), col = "red", from = .9, to = 2.1)
#' with(out$evals, points(x, fx))
#' for(k in 2:out$n_evals) {
#'   with(out$evals, abline(
#'     a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
#'     b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
#'     col = "blue"
#'   ))
#'   Sys.sleep(.5)
#' }
#'
#'
#'
#'
#'
#' f <- sin
#' x0 <- 2; x1 <- 1.9
#' (out <- unisecant(f, x0, x1))
#'
#' curve(f(x), col = "red", from = 0, to = 2*pi)
#' with(out$evals, points(x, fx))
#' for(k in 2:out$n_evals) {
#'   with(out$evals, abline(
#'     a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
#'     b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
#'     col = "blue"
#'   ))
#'   Sys.sleep(.5)
#' }
#'
#'
#'
#'
#' f <- log
#' x0 <- .40; x1 <- .50
#' (out <- unisecant(f, x0, x1))
#'
#' curve(f(x), col = "red", from = .25, to = 1.5)
#' with(out$evals, points(x, fx))
#' for(k in 2:out$n_evals) {
#'   with(out$evals, abline(
#'     a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
#'     b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
#'     col = "blue"
#'   ))
#'   Sys.sleep(.5)
#' }
#'
#'
#'
#' f <- function(x) (x-.24) * (x - .51) * (x - .76)
#' x0 <- .40; x1 <- .45
#' (out <- unisecant(f, x0, x1))
#'
#' curve(f(x), col = "red", from = .15, to = .85)
#' with(out$evals, points(x, fx))
#' for(k in 2:out$n_evals) {
#'   with(out$evals, abline(
#'     a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
#'     b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
#'     col = "blue"
#'   ))
#'   Sys.sleep(.5)
#' }
#'
#'
#'
#'
#' f <- function(x) pbeta(x, 6, 4) - .5
#' x0 <- .30; x1 <- .40
#' (out <- unisecant(f, x0, x1))
#'
#' curve(f(x), col = "red", from = 0, to = 1)
#' with(out$evals, points(x, fx))
#' for(k in 2:out$n_evals) {
#'   with(out$evals, abline(
#'     a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
#'     b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
#'     col = "blue"
#'   ))
#'   Sys.sleep(.5)
#' }
#'
#'
unisecant <- function(f, x0, x1, tol = 10*.Machine$double.eps, maxit = 100L) {

  # initialize x and fx
  x <- fx <- numeric(maxit)

  # check endpoint and return early if root there
  x[1]   <- x0;
  fx[1]  <- f(x[1])
  if(abs(fx[1]) <= tol) {
    return(list(
      root = x[1], val = fx[1],
      evals = data_frame(x = x[1], fx = fx[1]),
      n_evals = 1
    ))
  }

  x[2]   <- x1;
  fx[2]  <- f(x[2])
  if(abs(fx[2]) <= tol) {
    return(list(
      root = x[2], val = fx[2],
      evals = data_frame(x = x[1:2], fx = fx[1:2]),
      n_evals = 2
    ))
  }

  # loop
  for(k in 3:maxit) {
    x[k]  <- x[k-1] - fx[k-1] / ((fx[k-1]-fx[k-2])/(x[k-1]-x[k-2]))
    fx[k] <- f(x[k])
    n_evals <- k
    if(abs(fx[k]) <= tol) break
  }

  # return
  list(
    root = x[n_evals], f.root = fx[n_evals],
    evals = data_frame(x = x[1:n_evals], fx = fx[1:n_evals]),
    n_evals = n_evals
  )
}
