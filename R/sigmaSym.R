sigmaSym <-
function(X,Y,L,dimension,n,nvec,sigma){
  Id <- diag(1,dimension)
  npairs <- choose(nvec,2)
  Xdiff <- array(0,dim=c(max(npairs),dimension,L))
  enums <- array(0,dim=c(dimension,dimension,max(npairs),L))
  denoms <- array(0,dim=c(max(npairs),L))
  Xtmp <- array(0,dim=c(max(nvec),dimension,L))
  
  for (b in 1:L) {
    h <- 0
    Xtmp[1:dim(unique(X[which(Y==b),]))[1],,b] <- unique(X[which(Y==b),])
    for (j in 1:(nvec[b]-1)) {
      for (k in (j+1):nvec[b]) {
        h <- h+1
        Xdiff[h,,b] <- Xtmp[j,,b]-Xtmp[k,,b]
        enums[,,h,b] <-  Xdiff[h,,b] %*% t(Xdiff[h,,b])
      }
    }
  }
 
  # Fixed point iteration
  z <- 0
  repeat {
    sigma.old <- sigma
    sigma.inv <- solve(sigma)
    sigma <- matrix(0,dimension,dimension)
    for(b in 1:L) {
      denoms[1:npairs[b],b] <- apply( (Xdiff[1:npairs[b],,b] %*%
           sigma.inv) * Xdiff[1:npairs[b],,b],1,sum)
      sigma.tmp <- apply(sweep(enums[,,1:npairs[b],b],3,
                 denoms[1:npairs[b],b],"/"), c(1,2), sum)
      sigma <- sigma + sigma.tmp / nvec[b]
    }
    sigma <- dimension * 2 / (n-L) * sigma
    z <- z+1
    if(z>100) stop("did not converge!")
    if((Matrix::norm(solve(sigma.old)%*%sigma-diag(1,dimension),type="F")
        < 10^-5 | (z>100))) break
  }
  return(sigma)
}

