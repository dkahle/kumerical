#' Univariate bisection method
#'
#' Find the root of a function between two points a and b using the
#' bisection method.
#'
#' @param f function
#' @param a lower bound
#' @param b upper bound
#' @param tol tolerance, defaults to 10*.Machine$double.eps
#' @param maxit maximum number of iterations
#' @param bisections number or bisection steps for fast unibisect
#' @return a list
#' @name unibisect
#' @examples
#'
#' f <- function(x) x^2 - 2
#' a <- 0; b <- 2
#' unibisect(f, a, b)
#' simple_unibisect(f, a, b)
#' fast_unibisect(f, a, b)
#'
#' (out <- unibisect(f, a, b))
#'
#' curve(f(x), col = "red", from = a, to = b)
#' with(out$evals, points(x, fx))
#' out$n_evals # = number of function calls
#' plot(1:out$n_evals, out$evals$fx)
#'
#'
#'
#'
#' f <- sin
#' a <- .1; b <- 2*pi - .2
#' (out <- unibisect(f, a, b))
#'
#' curve(f(x), col = "red", from = a, to = b)
#' with(out$evals, points(x, fx))
#' out$n_evals # = number of function calls
#' plot(1:out$n_evals, out$evals$fx)
#'
#'
#'
#'
#' f <- function(x) (x-.24) * (x - .51) * (x - .76)
#' a <- 0; b <- 1
#' (out <- unibisect(f, a, b))
#'
#' curve(f(x), col = "red", from = a, to = b)
#' with(out$evals, points(x, fx))
#' out$n_evals # = number of function calls
#' plot(1:out$n_evals, out$evals$fx)
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




#' @export
#' @rdname unibisect
unibisect <- function(f, a, b, tol = 10*.Machine$double.eps, maxit = 100L) {

  # initialize x and fx
  x <- fx <- numeric(maxit)

  # check endpoints and return early if root there
  x[1] <- a;
  fa <- fx[1] <- f(a)
  if(abs(fa) <= tol) {
    return(list(
      root = a, f.root = fa,
      evals = data_frame(x = x[1], fx = fx[1]),
      n_evals = 1
    ))
  }

  x[2] <- b
  fb <- fx[2] <- f(b)
  if(sign(fa) == sign(fb)) stop("f(a) and f(b) must have opposite signs.")
  if(abs(fb) <= tol) {
    return(list(
      root = b, f.root = fb,
      evals = data_frame(x = x[1:2], fx = fx[1:2]),
      n_evals = 2
    ))
  }

  # compute first midpoint, return early if root found
  c  <- x[3]  <- (a+b)/2
  fc <- fx[3] <- f(c)
  if(abs(fc) < tol) {
    return(list(
      root = c, f.root = fc,
      evals = data_frame(x = x[1:3], fx = fx[1:3]),
      n_evals = 3
    ))
  }

  # loop
  for(k in 4:maxit) {
    if (sign(fa) == sign(fc)) {
      a <- c; fa <- fc
    } else {
      b <- c; fb <- fc
    }
    c  <- x[k]  <- (a+b)/2
    fc <- fx[k] <- f(c)
    n_evals <- k
    if(abs(fc) <= tol) break
  }

  # return
  if(abs(fc) > tol) warning("tolerance not achieved.")
  list(
    root = c, f.root = fc,
    evals = data_frame(x = x[1:n_evals], fx = fx[1:n_evals]),
    n_evals = n_evals
  )
}










#' @export
#' @rdname unibisect
fast_unibisect <- function(f, a, b, tol = 10*.Machine$double.eps, bisections = 8L, maxit = 100L) {

  if(bisections >= maxit) return(unibisect(f, a, b, tol, maxit))

  # initialize x and fx
  x <- fx <- numeric(maxit)

  # check endpoints and return early if root there
  x[1] <- a;
  fa <- fx[1] <- f(a)
  if(abs(fa) <= tol) {
    return(list(
      root = a, f.root = fa,
      evals = data_frame(x = x[1], fx = fx[1], method = "bisection"),
      n_evals = 1
    ))
  }

  x[2] <- b
  fb <- fx[2] <- f(b)
  if(sign(fa) == sign(fb)) stop("f(a) and f(b) must have opposite signs.")
  if(abs(fb) <= tol) {
    return(list(
      root = b, f.root = fb,
      evals = data_frame(x = x[1:2], fx = fx[1:2], method = "bisection"),
      n_evals = 2
    ))
  }

  # compute first midpoint, return early if root found
  c  <- x[3] <- (a+b)/2
  fc <- fx[3] <- f(c)
  if(abs(fc) < tol) {
    return(list(
      root = c, f.root = fc,
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
      root = c, f.root = fc,
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
  if(abs(fx[k]) > tol) warning("tolerance not achieved.")
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























#' @export
#' @rdname unibisect
simple_unibisect <- function(f, a, b, tol = 10*.Machine$double.eps, maxit = 100L) {

  # check endpoints and return early if root there
  fa <- f(a)
  if(abs(fa) <= tol) return(list(root = a, f.root = fa))

  fb <- f(b)
  if(sign(fa) == sign(fb)) stop("f(a) and f(b) must have opposite signs.")
  if(abs(fb) <= tol) return(list(root = b, f.root = fb))

  # compute first midpoint, return early if root found
  c  <- (a+b)/2; fc <- f(c)
  if(abs(fc) < tol) return(list(root = c, f.root = fc))

  # loop
  for(k in 4:maxit) {
    if (sign(fa) == sign(fc)) {
      a <- c; fa <- fc
    } else {
      b <- c; fb <- fc
    }
    c  <- (a+b)/2; fc <- f(c)
    if(abs(fc) <= tol) break
  }

  # return
  if(abs(fc) > tol) warning("tolerance not achieved.")
  list(root = c, f.root = fc)
}
