cvopttau <-
function(X,Y,th,Tau,wtype=c('linear','exponential'),distance=c('euclidean','ddeuclidean','mahalanobis'),cova=c('standard','M','sym'),nvec,n,L,dimension) {
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


  WX <- matrix(0,n,n)
  thIndex <- which(Y==th)
  HH <- rep(0,n)
  HH[thIndex] <- 1
  
  pi.old <- Inf
  
  for(tau in Tau){
    # compute T_th(X_i) for all i
    # Computation of the weight function W
    W <- rep(0,n)
    if(wtype=="exponential") {
      for(i in 1:n) {
        W[i] = (1-i/n)^tau
      }
    } else {
      if(wtype=="linear") {
        for(i in 1:n) {
          W[i] = pmax(1-(i/n)/tau,0)
        }
      }
    }
    for(i in 1:n){
      WX[i,order(di[i,])] <- W
    }
    # Compute test statistic
    for(i in 1:n) {
      T <- WX %*% HH
    }
    pi <- 0
    for(j in 1:n){
      #compute and sum pi_th(j)
      pi <- pi + sum(T[thIndex] <= T[j])

    }
    # min tau -> opt.tau
    if(pi < pi.old){
      pi.old <- pi
      opt.tau <- tau
    }
  }
  return(opt.tau)
}

  


