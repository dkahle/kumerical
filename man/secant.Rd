% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unisecant.R
\name{secant}
\alias{secant}
\alias{unisecant}
\alias{simple_unisecant}
\title{Univariate secant method}
\usage{
unisecant(f, x0, x1, tol = 10 * .Machine$double.eps, maxit = 100L)

simple_unisecant(f, x0, x1, tol = 10 * .Machine$double.eps, maxit = 100L)
}
\arguments{
\item{f}{function}

\item{x0}{initial value}

\item{x1}{second initial value}

\item{tol}{tolerance, defaults to 10*.Machine$double.eps}

\item{maxit}{maximum number of iterations}
}
\value{
a list
}
\description{
Find the root of a function using the secant method.
}
\examples{

f  <- function(x) x^2 - 2
x0 <- 0; x1 <- 1
unisecant(f, x0, x1)
simple_unisecant(f, x0, x1)

(out <- unisecant(f, x0, x1))

curve(f(x), col = "red", from = .9, to = 2.1)
with(out$evals, points(x, fx))
for(k in 2:out$n_evals) {
  with(out$evals, abline(
    a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
    b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
    col = "blue"
  ))
  Sys.sleep(.5)
}





f <- sin
x0 <- 2; x1 <- 1.9
(out <- unisecant(f, x0, x1))

curve(f(x), col = "red", from = 0, to = 2*pi)
with(out$evals, points(x, fx))
for(k in 2:out$n_evals) {
  with(out$evals, abline(
    a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
    b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
    col = "blue"
  ))
  Sys.sleep(.5)
}




f <- log
x0 <- .40; x1 <- .50
(out <- unisecant(f, x0, x1))

curve(f(x), col = "red", from = .25, to = 1.5)
with(out$evals, points(x, fx))
for(k in 2:out$n_evals) {
  with(out$evals, abline(
    a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
    b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
    col = "blue"
  ))
  Sys.sleep(.5)
}



f <- function(x) (x-.24) * (x - .51) * (x - .76)
x0 <- .40; x1 <- .45
(out <- unisecant(f, x0, x1))

curve(f(x), col = "red", from = .15, to = .85)
with(out$evals, points(x, fx))
for(k in 2:out$n_evals) {
  with(out$evals, abline(
    a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
    b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
    col = "blue"
  ))
  Sys.sleep(.5)
}




f <- function(x) pbeta(x, 6, 4) - .5
x0 <- .30; x1 <- .40
(out <- unisecant(f, x0, x1))

curve(f(x), col = "red", from = 0, to = 1)
with(out$evals, points(x, fx))
for(k in 2:out$n_evals) {
  with(out$evals, abline(
    a = fx[k-1] - (fx[k]-fx[k-1]) / (x[k]-x[k-1]) * x[k-1],
    b = (fx[k]-fx[k-1]) / (x[k]-x[k-1]),
    col = "blue"
  ))
  Sys.sleep(.5)
}

}
