pvs.gaussian <-
function(NewX,X,Y,cova=c('standard','M','sym')){
  # Adjust input
  X <- as.matrix(X)
  n<-dim(X)[1]
  dimension<-dim(X)[2]
  
  Y <- factor(Y)
  Y <- unclass(Y)
  L=max(Y)

  if(is.vector(NewX)) {
    NewX <- t(NewX)
  } else {
    NewX <- as.matrix(NewX)
  }
  
  nr<-dim(NewX)[1]
  s<-dim(NewX)[2]

  # Stop if dimensions of NewX[i,] and X[j,] do not match
  if(s!=dimension) {
    stop('dimensions of NewX[i,] and X[j,] do not match!')
  }

  PV <- matrix(rep(0,nr*L),nr)
  
  # Computation of nvec = vector with the numbers of training
  # observations from each class:
  nvec <- rep(0,L)
  for(b in 1:L) {
    nvec[b] <- sum(Y==b)
    
    # Stop if there are less than two observations from class b
    if(nvec[b]==0) {
      stop(paste('no observations from class',as.character(b),'!'))
    }
    if(nvec[b]==1) {
      stop(paste('only one observation from class', as.character(b),
                 '!'))
    }
  }
  
  # Compute mu
  mu <- matrix(0,L,dimension)
  for(m in 1:L) {
    mu[m,] = apply(X[Y==m,],2,mean)
  }
  cova <- match.arg(cova)
  
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
  
  for(i in 1:nr) {		
    #Add new observation NewX[i,] to Xtmp
    X.tmp <- rbind(X,NewX[i,])

    for(th in 1:L) {
      # Add NewX[i] temporarily to group th:
      thIndex <- c(which(Y==th),n+1)
      Y.tmp <- c(Y,th)
      nvec.tmp <- nvec
      nvec.tmp[th] <- nvec.tmp[th] + 1

      # Adjust mu
      mu.tmp <- mu
      mu.tmp[th,] <- mu[th,] + (X.tmp[n+1,]-mu[th,]) / (nvec[th]+1)
      
      #Adjust sigma
      if(cova=='standard'){
        sigma.tmp <- ((n-L)*sigma + 1/(1+1/nvec[th]) *
                     (NewX[i,]-mu[th]) %*%
                     t(NewX[i,]-mu[th])) / (n+1-L)
      } else {
        if(cova=='M'){
          sigma.tmp <- pvclass:::sigmaM(X=X.tmp,Y=Y.tmp,L=L,dimension=dimension,n=n+1,nvec=nvec.tmp,mu=mu.tmp,sigma=sigma)
        } else {
          if(cova=='sym'){
            sigma.tmp <- pvclass:::sigmaSym(X=X.tmp,Y=Y.tmp,L=L,dimension=dimension,n=n+1,nvec=nvec.tmp,sigma=sigma)
          } else {
            stop("type of covariance estimator not valid!")
          }
        }
      }
      sigma.tmp.inv <- solve(sigma.tmp)
      
      T <- matrix(0,n+1,1)
      for(j in thIndex){
        for(b in 1:L){
          T[j] <- T[j] + (b!=th)*nvec.tmp[b]/(n-nvec.tmp[th]) *
            exp((X.tmp[j,]-(mu.tmp[th,]+mu.tmp[b,])/2)%*%
                sigma.tmp.inv %*%(mu.tmp[b,]-mu.tmp[th,]))
        }
      }
      Tv <- T[thIndex]
     
      PV[i,th] <- sum(Tv >= Tv[length(Tv)]) / nvec.tmp[th]
    }
  }
  dimnames(PV)[[2]] <- attr(Y,'levels')
  return(PV)  
}

