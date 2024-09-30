hyperbootstrap <-
  function( xvec, B = 1000L, method = c("moments", "likelihood", "chisq"))
    
  {
    if (missing(xvec) )
      stop(paste("The arguments 'xvec' is missing.",
                 "Please, provide a vector of data in input."))
    nmm <- tolower(method)
    substr(nmm, 1, 1) <- toupper(substr(nmm, 1, 1))
    method <- match(method, c("moments", "likelihood", "chisq"))
    method <- na.omit(method)
    B <- as.integer(B)
    matsize <- B * length(method)
    res <- .C("hyperbootstrap", xvec, length(xvec), 
       a = double(matsize), b = double(matsize), 
       B, method, length(method), PACKAGE = "plpoisson")[3L:4L]
    res$a <- matrix(res$a, B, length(method))
    res$b <- matrix(res$b, B, length(method))
    colnames(res$a) <- colnames(res$b) <- nmm
   return(res)
  }
