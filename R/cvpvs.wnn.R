cvpvs.wnn <-
function(X,Y,wtype=c('linear','exponential'),W=NA,tau=0.3,distance=c('euclidean','ddeuclidean','mahalanobis'),cova=c('standard','M','sym')) {
  # Adjust input
  X <- as.matrix(X)
  n<-dim(X)[1]
  dimension<-dim(X)[2]

  Y <- factor(Y)
  Y <- unclass(Y)
  L=max(Y)
  
  wtype <- match.arg(wtype)
  cova <- match.arg(cova)
  distance <- match.arg(distance)

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
  
  # If tau is not a single value, search for the optimal tau.
  if( is.na(W)[1] & (length(tau)>1 | is.na(tau)[1]) ){
    if(wtype=='exponential'){
      # Use exponential weights
      if(is.na(tau)[1]) {
        tau <- c(1,5,10,20,30,40,50)
      }
    } else {
      if(wtype=='linear') {
        # Use linear weights
        if(is.na(tau)[1]) {
          tau <- seq(0.1,0.9,0.1)
        }
      }
    }
    
    opt.tau <- matrix(0,n,L)
    PV <- matrix(0,n,L)
    
    # optimal tau for correct classes
    for(th in 1:L){
      tmp <- pvclass:::cvopttau(X=X,Y=Y,th=th,Tau=tau,wtype=wtype,distance=distance,cova=cova,nvec=nvec,n=n,L=L,dimension=dimension)
      opt.tau[Y==th,th] <- tmp
      PV[Y==th,th] <- cvpvs.wnn(X=X,Y=Y,wtype=wtype,tau=tmp,distance=distance,cova=cova)[Y==th,th]
    }

    # optimal tau other classes
    for(i in 1:n){
      for(th in 1:L){
        if(Y[i]!=th){
          Y.tmp <- Y
          Y.tmp[i] <- th
          nvec.tmp <- nvec
          nvec.tmp[th] <- nvec[th]+1
          nvec.tmp[Y[i]] <- nvec[Y[i]]-1
          opt.tau[i,th] <- cvopttau(X=X,Y=Y.tmp,th=th,Tau=tau,wtype=wtype,distance=distance,cova=cova,nvec=nvec.tmp,n=n,L=L,dimension=dimension)
          PV[i,th] <- cvpvs.wnn(X=X,Y=Y,wtype=wtype,tau=opt.tau[i,th],distance=distance,cova=cova)[i,th]
        }
      }
    }
    attributes(PV)$opt.tau<-opt.tau
    dimnames(PV)[[2]] <- attr(Y,'levels')
    return(PV)
  }
  
 
  # Computation of the weight function W
  if(!is.na(W)[1]) {
    if(length(W)!=n) {
      stop(paste('length(W) != length(Y)'))
    }
  } else { 
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
  }
  
    
  WX <- matrix(0,n,n)
  
  if(distance=="mahalanobis") {
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
    
    for(i in 1:n) {
      Distance <- mahalanobis(X,X[i,],sigma)
      WX[i,order(Distance)] <- W
    }
    # Computation of the cross-validated P-values:
    PV <- matrix(0,n,L)

    # P-values for the correct classes, i.e. PV(i,Y[i]):
    AA <- matrix(0,n,1)
    for(th in 1:L) {
      thIndex <- which(Y==th)
      HH <- rep(0,n)
      HH[thIndex] <- 1
      AA <- WX %*% HH
      
      Tv <- AA[thIndex]
      PV[thIndex,th] <- pvclass:::CompPVs(Tv)
    }

    # P-values for other classes, i.e. PV(i,th), th ~= Y(i):
    for(i in 1:n) {
      for(th in 1:L) {
        if(th !=Y[i]) {
          
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
              sigma.tmp <- pvclass:::sigmaM(X=X,Y=Y.tmp,L=L,dimension=dimension,n=n,nvec=nvec.tmp,sigma=sigma)
            } else {
              if(cova=='sym'){
                sigma.tmp <- pvclass:::sigmaSym(X=X,Y=Y.tmp,L=L,dimension=dimension,n=n,nvec=nvec.tmp,sigma)
              } 
            }
          }
          
          thIndex <- c(i,which(Y==th))
          WX <- matrix(0,n,n)
          for(j in thIndex) {
            Distance <- mahalanobis(X,X[j,],sigma.tmp)
            WX[j,order(Distance)] <- W
          }
          
          HH <- rep(0,n)
          HH[thIndex] <- 1
          AA <- (WX %*% HH)
          
          Tv <- AA[thIndex]
          PV[i,th] <- sum(Tv <= Tv[1]) / nvec.tmp[th]
        }
      }		
    }

    dimnames(PV)[[2]] <- attr(Y,'levels')
    return(PV)
    
  } else {
    if(distance=="euclidean") {
      # Use Euclidean distance
      Distance <- as.matrix(dist(X))
      for(i in 1:n) {
        WX[i,order(Distance[i,])] <- W
      }
    } else {
      if(distance=="ddeuclidean") {
        # Use data driven Euclidean distance
        for(i in 1:dimension) {
          X[,i] <- X[,i]/sqrt(var(X[,i]))
        }
        
        Distance <- as.matrix(dist(X))
        for(i in 1:n) {
          WX[i,order(Distance[i,])] <- W
        }
      } else { 
        stop(paste("no valid distance!"))
      }
    }
  }

  # Computation of the cross-validated P-values:
  PV <- matrix(0,n,L)
  
  # P-values for the correct classes, i.e. PV(i,Y[i]):
  AA <- matrix(0,n,1)
  for(th in 1:L) {
    thIndex <- which(Y==th)
    HH <- rep(0,n)
    HH[thIndex] <- 1
    AA <- WX %*% HH
    
    Tv <- AA[thIndex]
    PV[thIndex,th] <- pvclass:::CompPVs(Tv)
  }
    
  # P-values for other classes, i.e. PV(i,th), th ~= Y(i):
  for(i in 1:n) {
    for(th in 1:L) {
      if(th !=Y[i]) {
        # Add observation i temporarily to class th
        thIndex <- c(i,which(Y==th))
        HH <- rep(0,n)
        HH[thIndex] <- 1
        AA <- (WX %*% HH)
          
        nvec_tmp <- nvec
        nvec_tmp[Y[i]] <- nvec[Y[i]] - 1
        nvec_tmp[th] <- nvec[th] + 1
        
        Tv <- AA[thIndex]
        PV[i,th] <- sum(Tv <= Tv[1]) / nvec_tmp[th]
      }
    }		
  }
  dimnames(PV)[[2]] <- attr(Y,'levels')
  return(PV)
}
