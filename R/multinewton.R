#' Univariate Newton method
#'
#' Find the root of a function using Newton's method.
#'
#' @param f function
#' @param df function; Jacobian matrix of f
#' @param x0 initial value
#' @param tol tolerance, defaults to 10*.Machine$double.eps
#' @param maxit maximum number of iterations
#' @return a list
#' @export
#' @examples
#'
#'
#' multinewton(function(x) x^2 - 2, function(x) 2*x, 2)
#'
#'
#'
#' # intersection of a plane, hyperboloid, and cone
#' f <- function(v) {
#'   x <- v[1]; y <- v[2]
#'   c(x^2 + y, x*y - 2*y^2)
#' }
#' df <- function(v) {
#'   x <- v[1]; y <- v[2]
#'   matrix(c(
#'     2*x,       1,
#'       y, x - 4*y
#'   ), nrow = 2, byrow = TRUE)
#' }
#' x0 <- c(2, 2)
#' f(x0)
#' df(x0)
#'
#' out <- multinewton(f, df, x0)
#' str(out, 1)
#'
#'
#'
#' # intersection of a plane, hyperboloid, and cone
#' # true solutions =
#' #   c(-3/sqrt(2), 0,  3/sqrt(2))
#' #   c( 3/sqrt(2), 0, -3/sqrt(2))
#' # corresponding to the nonlinear system
#' # x + y + z = 0
#' # x^2 - y^2 + z^2 = 9,
#' # x^2 + y^2 - z^2 = 0
#' f <- function(v) {
#'   x <- v[1]; y <- v[2]; z <- v[3]
#'   c(x + y + z, x^2 - y^2 + z^2 - 9, x^2 + y^2 - z^2)
#' }
#' df <- function(v) {
#'   x <- v[1]; y <- v[2]; z <- v[3]
#'   matrix(c(
#'       1,    1,    1,
#'     2*x, -2*y,  2*z,
#'     2*x,  2*y, -2*z
#'   ), nrow = 3, byrow = TRUE)
#' }
#' x0 <- c(2, 2, 2)
#' f(x0)
#' df(x0)
#' out <- multinewton(f, df, x0)
#' str(out, 1)
#' c( 3/sqrt(2), 0, -3/sqrt(2))
#' out$root
#'
#'
multinewton <- function(f, df, x0, tol = 10*.Machine$double.eps, maxit = 100L) {

  # initialize x and fx
  px <- length(x0)
  fx0 <- f(x0)
  pf <- length(fx0)
  x  <- matrix(NA_real_, nrow = maxit, ncol = px)
  fx <- matrix(NA_real_, nrow = maxit, ncol = pf)
  dfx <- replicate(maxit, matrix(NA_real_, nrow = pf, ncol = px), simplify = FALSE)

  # set norm
  norm <- function(v) sum(abs(v))

  # check endpoint and return early if root there
  x[1,]  <- x0;
  fx[1,] <- fx0
  if(norm(fx[1,]) <= tol) {
    return(list(
      root = x[1,], f.root = fx[1,],
      evals = data_frame(x = x[1,], fx = fx[1,], dfx = dfx[1]),
      n_evals = 1
    ))
  }

  # loop
  for(k in 2:maxit) {
    dfx[[k-1]] <- df(x[k-1,])
    h <- solve(dfx[[k-1]], -fx[k-1,])
    x[k,]  <- x[k-1,] + h
    fx[k,] <- f(x[k,])
    n_evals <- k
    if(norm(fx[k,]) <= tol) break
  }

  # return
  list(
    root = x[n_evals,], f.root = fx[n_evals,],
    evals = list("x" = x[1:n_evals,], "fx" = fx[1:n_evals,], "dfx" = dfx[1:(n_evals-1)]),
    n_evals = n_evals
  )
}
