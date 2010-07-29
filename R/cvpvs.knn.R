cvpvs.knn <-
function(X,Y,k=NA,distance=c('euclidean','ddeuclidean','mahalanobis'),cova=c('standard','M','sym')) {
  cova <- match.arg(cova)
  distance <- match.arg(distance)
  
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
  
  # If k is not a single integer, search for the optimal k.
  if(length(k)>1 | is.na(k)[1]) {
    if(is.na(k[1])) {
      k <- 1:ceiling(length(Y)/2)
    }
    opt.k <- matrix(0,n,L)
    PV <- matrix(0,n,L)
    
    # optimal k for correct classes
    for(th in 1:L){
      tmp <- pvclass:::cvoptk(X=X,Y=Y,th=th,K=k,distance=distance,cova=cova,nvec=nvec,n=n,L=L,dimension=dimension)
      opt.k[Y==th,th] <- tmp
      PV[Y==th,th] <- cvpvs.knn(X=X,Y=Y,k=tmp,distance=distance,cova=cova)[Y==th,th]
    }

    # optimal k other classes
    for(i in 1:n){
      for(th in 1:L){
        if(Y[i]!=th){
          Y.tmp <- Y
          Y.tmp[i] <- th
          nvec.tmp <- nvec
          nvec.tmp[th] <- nvec[th]+1
          nvec.tmp[Y[i]] <- nvec[Y[i]]-1
          opt.k[i,th] <- pvclass:::cvoptk(X=X,Y=Y.tmp,th=th,K=k,distance=distance,cova=cova,nvec=nvec.tmp,n=n,L=L,dimension=dimension)
          PV[i,th] <- cvpvs.knn(X=X,Y=Y,k=opt.k[i,th],distance=distance,cova=cova)[i,th]
        }
      }
    }
    attributes(PV)$opt.k<-opt.k
    dimnames(PV)[[2]] <- attr(Y,'levels')
    return(PV)
  }
  
   
  
  if(k>=n) stop('k >= sample size!')
  
  # Compute distances
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
    
    di <- matrix(0,n,n)
    for(i in 1:n) {
      di[i,] <- mahalanobis(X,X[i,],sigma)
    }
    # Computation of Rk :  a column vector, where Rk[i] is the
    #                      distance between X[i,] and its k-th nearest
    #                      neighbor.
    #                Nk :  Nk[i,b] is the number of training
    #                      observations from class b among the k
    #                      nearest neighbors of X[i,]
    Rk <- rep(0,n)
    Nk <- matrix(0,n,L)
    
    sdi <- t(apply(di,1,sort))
    Rk <- sdi[,k]
    for(i in 1:n) {
      for(b in 1:L) {
        bIndex <- which(Y==b)
        # get number of observations for each class among k
        # neighbours:
        Nk[i,b] <- sum(di[i,bIndex] <= Rk[i])
      }
    }
    
    # Computation of the cross-validated P-values:
    PV <- matrix(0,n,L)
    
    # P-values for the correct classes, i.e. PV[i,Y[i]]:
    for(th in 1:L) {
      JJ <- which(Y==th)
      Tv <- Nk[JJ,th]
      PV[JJ,th] <- pvclass:::CompPVs(Tv)
    }
    
    
    # P-values for other classes, i.e. PV[i,th], th != Y[i]:

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
          
          JJ <- c(i,which(Y==th))
          di <- matrix(0,n,n)
          for(j in JJ) {
            di[j,] <- mahalanobis(X,X[j,],sigma.tmp)
          }
          # Computation of Rk :  a column vector, where Rk[i] is the
          #                      distance between X[i,] and its k-th nearest
          #                      neighbor.
          #                Nk :  Nk[i,b] is the number of training
          #                      observations from class b among the k
          #                      nearest neighbors of X[i,]
          Rk <- rep(0,n)
          Nk <- matrix(0,n,L)
          
          sdi <- t(apply(di,1,sort))
          Rk <- sdi[,k]
          for(j in JJ) {
            Nk[j,th] <- sum(di[j,JJ] <= Rk[j])
          }
          
          Tv <- Nk[JJ,th]
          PV[i,th] <- sum(Tv <= Tv[1]) / nvec.tmp[th]
        }
      }
    }
    dimnames(PV)[[2]] <- attr(Y,'levels')
    return(PV)
  } else {
    # Use Euclidean distance

    if(distance=="ddeuclidean") {
      # Use data driven Euclidean distance
      for(i in 1:dimension) {
        X[,i] <- X[,i]/sqrt(var(X[,i]))
      }
    }
    di  <- as.matrix(dist(X))
  }
  
  # Computation of Rk :  a column vector, where Rk[i] is the
  #                      distance between X[i,] and its k-th nearest
  #                      neighbor.
  #                Nk :  Nk[i,b] is the number of training
  #                      observations from class b among the k
  #                      nearest neighbors of X[i,]
  Rk <- rep(0,n)
  Nk <- matrix(0,n,L)
  
  sdi <- t(apply(di,1,sort))
  Rk <- sdi[,k]
  for(i in 1:n) {
    for(b in 1:L) {
      bIndex <- which(Y==b)
      # get number of observations for each class among k
      # neighbours:
      Nk[i,b] <- sum(di[i,bIndex] <= Rk[i])
    }
  }
  
  # Computation of the cross-validated P-values:
  PV <- matrix(0,n,L)
    
  # P-values for the correct classes, i.e. PV[i,Y[i]]:
  for(th in 1:L) {
    JJ <- which(Y==th)
    Tv <- Nk[JJ,th]
    PV[JJ,th] <- pvclass:::CompPVs(Tv)
  }
    
    
  # P-values for other classes, i.e. PV[i,th], th ~= Y[i]:
  for(i in 1:n) {
    for(th in 1:L) {
      if(th !=Y[i]) {
        JJ <- c(i,which(Y==th))
        Nk_tmp  <- Nk
        nvec_tmp <- nvec
        Nk_tmp[JJ,th] <- Nk_tmp[JJ,th] + (di[JJ,i] <= Rk[JJ])
        Nk_tmp[JJ,Y[i]] <- Nk_tmp[JJ,Y[i]] - (di[JJ,i] <= Rk[JJ])
        nvec_tmp[Y[i]] <- nvec[Y[i]] - 1
        nvec_tmp[th] <- nvec[th] + 1
        
        Tv <- Nk_tmp[JJ,th]
        PV[i,th] <- sum(Tv <= Tv[1]) / nvec_tmp[th]
      }
    }
  }

  dimnames(PV)[[2]] <- attr(Y,'levels')
  return(PV)
}

