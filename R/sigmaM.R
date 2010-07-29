sigmaM <-
function(X,Y,L,dimension,n,nvec,mu=NA,sigma=NA){
  if(any(is.na(mu))){
    # Compute mu
    mu <- matrix(0,L,dimension)
    for(m in 1:L) {
      mu[m,] = apply(X[Y==m,],2,mean)
    }
  }

  M <- array(0,c(dimension,dimension,L))
	Xc <- X - mu[Y,]
	for (b in 1:L){
		M[,,b] <- t(Xc[Y==b,]) %*% Xc[Y==b,]
		}

  if(any(is.na(sigma))){
    # Compute sigma
    sigma <- apply(M,c(1,2),sum)/(n-L)
  }
  # Fixed point iteration
  repeat{
    sigma.old <- sigma
    sigma <- matrix(0,dimension,dimension)
    for(b in 1:L){
      sigma <- sigma + nvec[b] * M[,,b] /
        sum(diag( solve(sigma.old) %*% M[,,b] ))
    }
    sigma <- dimension / n * sigma
    if(Matrix::norm(solve(sigma.old)%*%sigma - diag(1,dimension),type="F")    
        < 10^-5) break
  }
  return(sigma)
}

