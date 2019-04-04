### Forbes index functions by J. Alroy ###
##### Available from http://bio.mq.edu.au/~jalroy/Forbes.R ###
##### help page http://bio.mq.edu.au/~jalroy/Forbes.html ###

forbes<-function(x,y,corrected)	{
  if (missing(corrected))	{
    corrected <- T
  }
  if (is.numeric(x) && is.numeric(y) && min(x) == 0 && min(y) == 0 && length(x) == length(y))	{
    a <- length(which((x * y) > 0))
    b <- length(which(x > 0)) - a
    c <- length(which(y > 0)) - a
  } else	{
    a <- length(na.omit(match(x,y)))
    b <- length(x) - a
    c <- length(y) - a
  }
  n <- a + b + c
  if (corrected == T)	{
    return(a * (n + sqrt(n))/((a + b) * (a + c) + a * sqrt(n) + (b * c)/2))
  } else	{
    return(a * n/((a + b) * (a + c)))
  }
}

forbesMatrix<-function(x,corrected)	{
  x[is.na(x)] <- 0
  m = matrix(nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x))	{
    for (j in 1:ncol(x))	{
      m[i,j] = forbes(x[,i],x[,j],corrected)
    }
  }
  m <- as.data.frame(m)
  rownames(m) <- colnames(x)
  colnames(m) <- colnames(x)
  return(m)
}
