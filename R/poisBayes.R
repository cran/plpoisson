poisBayes <-
function( xobs, n, s, t, a, b, alpha = 0.05)

{
  if (missing(a) || missing(b))
    stop(paste("The arguments 'a' and 'b' are missing.",
         "Please, provide values in input."))
  ns <- n * s
  nst <- ns + t
  if (is.finite(b)) {
    prob <- (b * ns + 1) / (b * nst + 1)
  } else { # when 'b' is infinite
    prob <- ns / nst
  }

  res <- list()
  res$lower <- qnbinom(alpha, xobs + a, prob)
  res$upper <- qnbinom(alpha, xobs + a, prob, lower.tail = FALSE)
  return(res)
}
