poisUNIF <-
function(xobs, n, s, t, alpha = 0.05)

{
  ns <- n * s
  nst <- ns + t
  prob <- ns / nst

  res <- list()
  res$lower <- qnbinom(alpha, xobs + 1, prob)
  res$upper <- qnbinom(alpha, xobs + 1, prob, lower.tail = FALSE)
  return(res)
}
