poisJEFF <-
function(xobs, n, s, t, alpha = 0.05)

{
  ns <- n * s
  nst <- ns + t
  prob <- ns / nst

  res <- list()
  res$lower <- qnbinom(alpha, xobs, prob)
  res$upper <- qnbinom(alpha, xobs, prob, lower.tail = FALSE)
  return(res)
}
