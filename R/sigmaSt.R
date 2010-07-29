sigmaSt <-
function(X,Y,L,dimension,n,mu=NA){
  if(any(is.na(mu))){
    # Compute mu
    mu <- matrix(0,L,dimension)
    for(m in 1:L) {
      mu[m,] = apply(X[Y==m,],2,mean)
    }
  }
  
  #Compute sigma
  Xc <- X - mu[Y,]
  sigma <- t(Xc) %*% Xc / (n-L)
}

