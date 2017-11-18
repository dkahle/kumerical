#' Univariate Newton method
#'
#' Find the root of a function using Newton's method.
#'
#' @param f function
#' @param df function; derivative of f
#' @param x0 initial value
#' @param tol tolerance, defaults to 10*.Machine$double.eps
#' @param maxit maximum number of iterations
#' @return a list
#' @export
#' @examples
#'
#'
#' uninewton(function(x) x^2 - 2, function(x) 2*x, 2)
#'
#'
#' f  <- function(x) x^2 - 2
#' df <- function(x) 2*x
#' x0 <- 2
#' (out <- uninewton(f, df, x0))
#'
#' curve(f(x), col = "red", from = .9, to = 2.1)
#' with(out$evals, points(x, fx))
#' for(k in 1:out$n_evals) {
#'   with(out$evals, abline(a = fx[k] - dfx[k]*x[k], b = dfx[k], col = "blue"))
#'   Sys.sleep(1)
#' }
#'
#'
#'
#'
#'
#' f <- sin
#' df <- cos
#' x0 <- 2
#' (out <- uninewton(f, df, x0))
#'
#' curve(f(x), col = "red", from = 0, to = 2*pi)
#' with(out$evals, points(x, fx))
#' for(k in 1:out$n_evals) {
#'   with(out$evals, abline(a = fx[k] - dfx[k]*x[k], b = dfx[k], col = "blue"))
#'   Sys.sleep(1)
#' }
#'
#'
#'
#'
#' f <- log
#' df <- function(x) 1/x
#' x0 <- .40
#' (out <- uninewton(f, df, x0))
#'
#' curve(f(x), col = "red", from = .25, to = 1.5)
#' with(out$evals, points(x, fx))
#' for(k in 1:out$n_evals) {
#'   with(out$evals, abline(a = fx[k] - dfx[k]*x[k], b = dfx[k], col = "blue"))
#'   Sys.sleep(1)
#' }
#'
#'
#'
#' f <- function(x) (x-.24) * (x - .51) * (x - .76)
#' df <- function(x) 3*x^2 - 3.02*x + .6924
#' x0 <- .40
#' (out <- uninewton(f, df, x0))
#'
#' curve(f(x), col = "red", from = .15, to = .85)
#' with(out$evals, points(x, fx))
#' for(k in 1:out$n_evals) {
#'   with(out$evals, abline(a = fx[k] - dfx[k]*x[k], b = dfx[k], col = "blue"))
#'   Sys.sleep(1)
#' }
#'
#'
uninewton <- function(f, df, x0, tol = 10*.Machine$double.eps, maxit = 100L) {

  # initialize x and fx
  x <- fx <- dfx <- numeric(maxit)

  # check endpoint and return early if root there
  x[1]   <- x0;
  fx[1]  <- f(x[1])
  if(abs(fx[1]) <= tol) {
    return(list(
      root = x[1], f.root = fx[1],
      evals = data_frame(x = x[1], fx = fx[1], dfx = NA),
      n_evals = 1
    ))
  }

  # loop
  for(k in 2:maxit) {
    dfx[k-1] <- df(x[k-1])
    x[k]  <- x[k-1] - fx[k-1] / dfx[k-1]
    fx[k] <- f(x[k])
    n_evals <- k
    if(abs(fx[k]) <= tol) break
  }

  # return
  list(
    root = x[n_evals], f.root = fx[n_evals],
    evals = data_frame("x" = x[1:n_evals], "fx" = fx[1:n_evals], "dfx" = dfx[1:n_evals]),
    n_evals = n_evals
  )
}
