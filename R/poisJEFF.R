poisJEFF <-
function(xobs, n, s, t, alpha = 0.05)

{
  ns <- n * s
  nst <- ns + t
  prob <- ns / nst

  res <- list()
  res$lower <- qnbinom(alpha, xobs, prob) + 1
  res$upper <- qnbinom(alpha, xobs, prob, lower.tail = FALSE) - 1
  return(res)
}
