cvoptk <-
function(X,Y,th,K,distance=c('euclidean','ddeuclidean','mahalanobis'),cova=c('standard','M','sym'),nvec,n,L,dimension) {
  
  if(distance=='mahalanobis'){
    # Use Mahalanobis distance
    # Compute mu
    mu <- matrix(0,L,dimension)
    for(m in 1:L) {
      mu[m,] = apply(X[Y==m,],2,mean)
    }
    
    # Compute sigma      
    if(cova=='standard'){
      sigma <- pvclass:::sigmaSt(X=X,Y=Y,L=L,dimension=dimension,n=n,mu=mu)
    } else {
      if(cova=='M'){
        sigma <- pvclass:::sigmaM(X=X,Y=Y,L=L,dimension=dimension,n=n,nvec=nvec,mu=mu)
      } else {
        if(cova=='sym'){
          sigma <- pvclass:::sigmaSt(X=X,Y=Y,L=L,dimension=dimension,n=n,mu=mu)
          sigma <- pvclass:::sigmaSym(X=X,Y=Y,L=L,dimension=dimension,n=n,nvec=nvec,sigma=sigma)
        } else {
          stop("type of covariance estimator not valid!")
        }
      }
    }
    
    di <- matrix(0,n,n)
    for(i in 1:n) {
      di[i,] <- mahalanobis(X,X[i,],sigma)
    }
  } else{
    if(distance=='ddeuclidean'){
      # Use data driven Euclidean distance
      for(i in 1:dimension) {
        X[,i] <- X[,i]/sqrt(var(X[,i]))
      }
    }
    di  <- as.matrix(dist(X))
  }

  
  Rk <- rep(0,n)
  Nk <- Rk
  sdi <- t(apply(di,1,sort))
  thIndex <- which(Y==th)
  pi.old <- Inf

  for(k in K){
    #compute T_th(X_i) for all i
    Rk <- sdi[,k]
    for(i in 1:n) {
      Nk[i] <- sum(di[i,thIndex] <= Rk[i])
    }
    pi <- 0
    for(j in 1:n){
      #compute and sum pi_th(j)
      pi <- pi + sum(Nk[thIndex] <= Nk[j])
    }
    # min k -> opt.k
    if(pi < pi.old){
      pi.old <- pi
      opt.k <- k
    }
  }
  return(opt.k)
}

  


