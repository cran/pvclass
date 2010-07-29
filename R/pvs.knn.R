pvs.knn <-
function(NewX,X,Y,k=NA,distance=c('euclidean','ddeuclidean','mahalanobis'),cova=c('standard','M','sym')){
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
  
  cova <- match.arg(cova)

  # If k is not a single integer, choose optimal k
  if(length(k)>1 | is.na(k)[1]) {
    if(is.na(k[1])) {
      k <- 1:ceiling(length(Y)/2)
    }
    opt.k <- matrix(0,nr,L)
    for(i in 1:nr){
      for(j in 1:L){
        sum.pi <- rep(0,length(k))
        for(l in 1:length(k)){
          a <- cvpvs.knn(X=rbind(X,NewX[i,]),Y=c(Y,j),k=k[l],distance=distance,cova=cova)
          b <- 0
          for(m in 1:L){
            b <- b + sum(a)
          }
          sum.pi[l] <- b
        }
        opt.k[i,j] <- k[which.min(sum.pi)]
      }
    }
    PV <- matrix(0,nr,L)
    for(i in 1:nr){
      for(j in 1:L){
        PV[i,j] <- pvs.knn(NewX=NewX[i,],X=X,Y=Y,k=opt.k[i,j],distance=distance,cova=cova)[j]
      }
    }
    dimnames(opt.k)[[2]] <- attr(Y,'levels')
    dimnames(PV)[[2]] <- attr(Y,'levels')
    attributes(PV)$opt.k<-opt.k
    return(PV)
  }
  
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

  
  distance <- match.arg(distance)
  if(distance=="mahalanobis") {
    # Use Mahalanobis distance
      
    di <- matrix(0,n+1,n+1)
    
    # Compute mu
    mu <- matrix(0,L,dimension)
    for(m in 1:L) {
      mu[m,] = apply(X[Y==m,],2,mean)
    }
    
    # Compute sigma      
    if(cova=='standard'){
      sigma <- sigmaSt(X=X,Y=Y,L=L,dimension=dimension,n=n,mu=mu)
    } else {
      if(cova=='M'){
        sigma <- sigmaM(X=X,Y=Y,L=L,dimension=dimension,n=n,nvec=nvec,mu=mu)
      } else {
        if(cova=='sym'){
          sigma <- sigmaSt(X=X,Y=Y,L=L,dimension=dimension,n=n,mu=mu)
          sigma <- sigmaSym(X=X,Y=Y,L=L,dimension=dimension,n=n,nvec=nvec,sigma=sigma)
        } else {
          stop("type of covariance estimator not valid!")
        }
      }
    }
    Nk <- matrix(0,n+1,1)
    
    for(i in 1:nr) {		
      # Add new observation NewX[i,] to Xtmp
      Xtmp <- rbind(X,NewX[i,])
      
      for(th in 1:L) {
        # Add NewX[i] temporarily to group th:
        thIndex <- c(which(Y==th),n+1)
        Ytmp <- c(Y,th)
        nvec.tmp <- nvec
        nvec.tmp[th] <- nvec.tmp[th] + 1
        
        # Adjust sigma
        if(cova=='standard'){
          sigmatmp <- ((n-L)*sigma + 1/(1+1/nvec[th]) *
                       (NewX[i,]-mu[th]) %*%
                       t(NewX[i,]-mu[th])) / (n+1-L)
        } else {
          if(cova=='M'){
            sigmatmp <- sigmaM(X=Xtmp,Y=Ytmp,L=L,dimension=dimension,n=n+1,nvec=nvec.tmp,sigma=sigma)
          } else {
            if(cova=='sym'){
              sigmatmp <- sigmaSym(X=Xtmp,Y=Ytmp,L=L,dimension=dimension,n=n+1,nvec=nvec.tmp,sigma=sigma)
            } else {
              stop("type of covariance estimator not valid!")
            }
          }
        }
        
        
        # Compute Mahalanobis distances
        for(m in 1:(n+1)) {
          di[m,] <- mahalanobis(Xtmp,Xtmp[m,],sigmatmp)
        }
        
        sdi <- t(apply(di,1,sort))
        Rk <- sdi[,k]
        for(j in 1:(n+1)) {
          Nk[j] <- sum(di[j,thIndex]<=Rk[j])            
        }
        

        # Compute PV[i,th]
        Tv <- Nk[thIndex]
        PV[i,th] <- sum(Tv <= Tv[length(Tv)]) / nvec.tmp[th]
        
      }
    }
    dimnames(PV)[[2]] <- attr(Y,'levels')
    return(PV)  
          
    
    
    } else {
      if(distance=="ddeuclidean") {
        # Use data driven Euclidean distance
        Nk <- matrix(0,n+1,1)
        for(i in 1:nr) {		
          # Add new observation NewX[i,] to Xtmp
          Xtmp <- rbind(X,NewX[i,])

          # Compute data driven euclidean distances
          for(m in 1:dimension) {
            Xtmp[,m] <- Xtmp[,m]/sqrt(var(Xtmp[,m]))
          }
          
          di <- as.matrix(dist(Xtmp))
          sdi <- t(apply(di,1,sort))
          Rk <- sdi[,k]
          
          for(th in 1:L) {
            # Add NewX[i] temporarily to group th:
            thIndex <- c(which(Y==th),n+1)
            nvec_tmp <- nvec
            nvec_tmp[th] <- nvec_tmp[th] + 1
            
            for(j in 1:(n+1)) {
              Nk[j] <- sum(di[j,thIndex]<=Rk[j])            
            }
          
            
            # Compute PV[i,th]
            Tv <- Nk[thIndex]
            PV[i,th] <- sum(Tv <= Tv[length(Tv)]) / nvec_tmp[th]
          }
        }
        dimnames(PV)[[2]] <- attr(Y,'levels')
        return(PV) 
      } else {
        # Use Euclidean distance
        di <- as.matrix(dist(rbind(X,NewX)))
        Nk <- matrix(0,n,L)
        Nkm1 <- matrix(0,n,L)
        sdi <- t(apply(di[1:n,1:n],1,sort))
        Rk <- sdi[,k]
        rkm1 <- sdi[,k-1]
        for(j in 1:n) {
          for(th in 1:L) {
            thIndex <- which(Y==th)
            Nk[j,th] <- sum(di[1:n,1:n][j,thIndex]<=Rk[j])
            Nkm1[j,th] <- sum(di[1:n,1:n][j,thIndex]<=rkm1[j])
          }
        }
        for(i in 1:nr) {		
          Nk_tmp <- rbind(cbind(Nk,rep(0,n)),rep(0,L+1))
          Nk_tmp[n+1,L+1] <- 1
          # Add new observation NewX[i,] to Xtmp
          Xtmp <- rbind(X,NewX[i,])
          di_tmp <- di[c(1:n,n+i),c(1:n,n+i)]
          sdi_tmp <- t(apply(di_tmp,1,sort))
          NewRk <- sdi_tmp[,k]
          
          # Determine Nk_tmp[n+1,1:L]
          for(th in 1:L) {
            thIndex <- which(Y==th)
            Nk_tmp[n+1,th] <- sum(di_tmp[n+1,thIndex]<=NewRk[n+1])
          }
          # Update Nk_tmp[1:n,]
          JJ1 <- which(di_tmp[n+1,1:n]<Rk)
          JJ2 <- which(di_tmp[n+1,1:n]==Rk)
          Nk_tmp[JJ1,1:L] <- Nkm1[JJ1,]
          Nk_tmp[JJ1,L+1] <- 1
          Nk_tmp[JJ2,L+1] <- 1
         # Nk_tmp <- matrix(0,n+1,L)
          for(th in 1:L) {
            # Add NewX[i] temporarily to group th:
            thIndex <- c(which(Y==th),n+1)
            Nk2 <- Nk_tmp[thIndex,1:L]
            Nk2[,th] <- Nk2[,th] + Nk_tmp[thIndex,L+1]
            nvec_tmp <- nvec
            nvec_tmp[th] <- nvec_tmp[th] + 1
            

            # Compute PV[i,th]
            Tv <- Nk2[,th]
            PV[i,th] <- sum(Tv <= Tv[length(Tv)]) / nvec_tmp[th]
          }
        }
        dimnames(PV)[[2]] <- attr(Y,'levels')
        return(PV)
      }
    }
}
