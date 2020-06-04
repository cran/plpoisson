poiss <-
function(xobs, n, s, t, alpha = 0.05)

{
  ns <- n * s
  nst <- ns + t
  prob <- ns / nst

  res <- .C("plBinom", as.integer(xobs), prob, alpha, 
            upper = -1L, lower = -1L, 
            PACKAGE = "plpoisson")[5L:4L]
  return(res)
}
