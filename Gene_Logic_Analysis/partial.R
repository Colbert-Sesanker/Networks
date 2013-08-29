# For partial application of functions, partial(function, arg1,...,argn) returns a function with supplied arguments bound, and an airity of the
# number of arguments not supplied
# Colbert Sesanker Jan 2013

partial <- function(f, ...) {
  capture <- list(...)
  function(...) {
    do.call(f, c(capture, list(...)))
  }
}
