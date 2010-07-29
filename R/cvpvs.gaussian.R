cvpvs.gaussian <-
function(X,Y,cova=c('standard','M','sym')) {
  cova <- match.arg(cova)

  # Adjust input
  X <- as.matrix(X)
  n<-dim(X)[1]
  dimension<-dim(X)[2]

  Y <- factor(Y)
  Y <- unclass(Y)
  L=max(Y)
 
  # Computation of nvec: an L-dimensional column vector, where nvec[b]
  # is the number of training observations belonging to class b.
  nvec <- rep(0,L)
  for(b in 1:L) {
    nvec[b] <- sum(Y==b)
    
    # Stop if there are less than two observations from class b
    if(nvec[b]==0) {
      stop(paste('no observations from class',
                 as.character(b),'!'))
    }
    if(nvec[b]==1) {
      stop(paste('only one observation from class',
                 as.character(b),'!'))
    }
  }
  
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
  sigma.inv <- solve(sigma)
  
  # Computation of the cross-validated P-values:
  PV <- matrix(0,n,L)
    
  # P-values for the correct classes, i.e. PV[i,Y[i]]:
    
  T <- matrix(0,n,L)
  for(i in 1:n){
    for(th in 1:L){
      for(b in 1:L){
        T[i,th] <- T[i,th] + nvec[b]/(n-nvec[th]) *
          exp(X[i,]-(mu[th,]+mu[b,])/2)%*% sigma.inv %*%
            (mu[b,]-mu[th,])
      }
    }
  }

  for(th in 1:L) {
    JJ <- which(Y==th)
    Tv <- -T[JJ,th]
    PV[JJ,th] <- pvclass:::CompPVs(Tv)
  }
  
    
  # P-values for other classes, i.e. PV[i,th], th ~= Y[i]:
  for(i in 1:n) {
    for(th in 1:L) {
      if(th !=Y[i]) {
        
        T <- matrix(0,n,L)
        # Adjust mu
        mu.tmp <- mu
        mu.tmp[Y[i],] <- mu[Y[i],] - (X[i,]-mu[Y[i],]) / (nvec[Y[i]]-1)
        mu.tmp[th,] <- mu[th,] + (X[i,]-mu[th,]) / (nvec[th]+1)

        # Adjust nvec
        nvec.tmp <- nvec
        nvec.tmp[c(Y[i],th)] <- c(nvec[Y[i]]-1,nvec[th]+1)
        
        # Adjust sigma
        if(cova=='standard'){
          sigma.tmp <- ((n-L)*sigma - (X[i,]-mu[Y[i],])%*%
                        t(X[i,]-mu[Y[i],]) / 1-1/nvec[Y[i]] +
                        (X[i,]-mu[th,])%*%t(X[i,]-mu[th,]) /
                        1-1/nvec[th]) / (n-L)
        } else {
          Y.tmp <- Y
          Y.tmp[i] <- th
          if(cova=='M'){
            sigma.tmp <- pvclass:::sigmaM(X=X,Y=Y.tmp,L=L,dimension=dimension,n=n,nvec=nvec.tmp,mu=mu.tmp,sigma=sigma)
          } else {
            if(cova=='sym'){
              sigma.tmp <- pvclass:::sigmaSym(X=X,Y=Y.tmp,L=L,dimension=dimension,n=n,nvec=nvec.tmp,sigma)
            } 
          }
        }
        sigma.tmp.inv <- solve(sigma.tmp)
                
        # Compute T
        JJ <- c(i,which(Y==th))
        for(j in JJ){
          for(b in 1:L){
            T[j,th] <- T[j,th] + (b!=th)*nvec.tmp[b]/(n-nvec.tmp[th]) *
              exp((X[j,]-(mu.tmp[th,]+mu.tmp[b,])/2)%*%
                sigma.tmp.inv %*%(mu.tmp[b,]-mu.tmp[th,]))
          }
        }
               
        Tv <- T[JJ,th]
        PV[i,th] <- sum(Tv >= Tv[1]) / nvec.tmp[th]
      }
    }
  }
  dimnames(PV)[[2]] <- attr(Y,'levels')
  return(PV)
}

