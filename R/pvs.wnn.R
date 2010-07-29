pvs.wnn <-
function(NewX,X,Y,wtype=c('linear','exponential'),W=NA,tau=0.3,distance=c('euclidean','ddeuclidean','mahalanobis'),cova=c('standard','M','sym')){
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
  wtype <- match.arg(wtype)

  # If tau is not a single value choose optimal tau
  if( is.na(W)[1] & (length(tau)>1 | is.na(tau)[1]) ){
    if(wtype=='exponential'){
      # Use exponential weights
      if(is.na(tau)[1]) {
        tau <- c(1,5,10,20,30,40,50)
      }
      opt.tau <- matrix(0,nr,L)
      for(i in 1:nr){
        for(j in 1:L){
          sum.pi <- rep(0,length(tau))
          for(l in 1:length(tau)){
            a <- cvpvs.wnn(X=rbind(X,NewX[i,]),Y=c(Y,j),wtype="e",tau=tau[l],distance=distance,cova=cova)
            b <- 0
            for(m in 1:L){
              b <- b + sum(a)
            }
            sum.pi[l] <- b
          }
          opt.tau[i,j] <- tau[which.min(sum.pi)]
        }
      }
      PV <- matrix(0,nr,L)
      for(i in 1:nr){
        for(j in 1:L){
          PV[i,j] <- pvs.wnn(NewX=NewX[i,],X=X,Y=Y,wtype=wtype,tau=opt.tau[i,j],distance=distance,cova=cova)[j]
        }
      }
      dimnames(opt.tau)[[2]] <- attr(Y,'levels')
      dimnames(PV)[[2]] <- attr(Y,'levels')
      attributes(PV)$opt.tau<-opt.tau
      return(PV)
      
    } else {
      if(wtype=='linear') {
        # Use linear weights
        if(is.na(tau)[1]) {
          tau <- seq(0.1,0.9,0.1)
        }
        opt.tau <- matrix(0,nr,L)
        for(i in 1:nr){
          for(j in 1:L){
            sum.pi <- rep(0,length(tau))
            for(l in 1:length(tau)){
              a <- cvpvs.wnn(X=rbind(X,NewX[i,]),Y=c(Y,j),wtype="l",tau=tau[l],distance=distance,cova=cova)
              b <- 0
              for(m in 1:L){
                b <- b + sum(a)
              }
              sum.pi[l] <- b
            }
            opt.tau[i,j] <- tau[which.min(sum.pi)]
          }
        }
        PV <- matrix(0,nr,L)
        for(i in 1:nr){
          for(j in 1:L){
            PV[i,j] <- pvs.wnn(NewX=NewX[i,],X=X,Y=Y,wtype=wtype,tau=opt.tau[i,j],distance=distance,cova=cova)[j]
          }
        }
        dimnames(opt.tau)[[2]] <- attr(Y,'levels')
        dimnames(PV)[[2]] <- attr(Y,'levels')
        attributes(PV)$opt.tau<-opt.tau
        return(PV)
      }
    }
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
  
  # Computation of the weight function W
  if(!is.na(W)[1]) {
    if(length(W)!=(n+1)) {
      stop(paste('length(W) != length(Y)+1'))
    }
  } else { 
    W <- rep(0,n+1)
    if(wtype=="exponential") {
      for(i in 1:(n+1)) {
        W[i] = (1-i/(n+1))^tau
      }
    } else {
      if(wtype=="linear") {
        for(i in 1:(n+1)) {
          W[i] = pmax(1-(i/(n+1))/tau,0)
        }
      }
    }
  }
  
  distance <- match.arg(distance)
  if(distance=="mahalanobis") {
    # Use Mahalanobis distance
    WXtmp <- matrix(0,n+1,n+1)
    di <- matrix(0,n+1,n+1)
    
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
    
    for(i in 1:nr) {		
      # Add new observation NewX[i,] to Xtmp
      Xtmp <- rbind(X,NewX[i,])
      
      for(th in 1:L) {
        # Add NewX[i] temporarily to group th:
        thIndex <- c(which(Y==th),n+1)
        Ytmp <- c(Y,th)
        HH <- rep(0,n+1)
        HH[thIndex] <- 1
        nvec.tmp <- nvec
        nvec.tmp[th] <- nvec.tmp[th] + 1
        
        # Adjust sigma
        if(cova=='standard'){
          sigmatmp <- ((n-L)*sigma + 1/(1+1/nvec[th]) *
                       (NewX[i,]-mu[th]) %*%
                       t(NewX[i,]-mu[th])) / (n+1-L)
        } else {
          if(cova=='M'){
            sigmatmp <- pvclass:::sigmaM(X=Xtmp,Y=Ytmp,L=L,dimension=dimension,n=n+1,nvec=nvec.tmp,sigma=sigma)
          } else {
            if(cova=='sym'){
              sigmatmp <- pvclass:::sigmaSym(X=Xtmp,Y=Ytmp,L=L,dimension=dimension,n=n+1,nvec=nvec.tmp,sigma=sigma)
            } else {
              stop("type of covariance estimator not valid!")
            }
          }
        }
        
        # Compute Mahalanobis distances
        for(m in 1:(n+1)) {
          di[m,] <- mahalanobis(Xtmp,Xtmp[m,],sigmatmp)
        }
        
        # Assign weights
        for(j in 1:(n+1)) {
          WXtmp[j,order(c(di[j,c(1:(n+1))]))] <- W
        }
                
        # Compute PV[i,th]
        AA <- (WXtmp %*% HH)
        Tv <- AA[thIndex]
        PV[i,th] <- sum(Tv <= Tv[length(Tv)]) / nvec.tmp[th]
      }
    }
    dimnames(PV)[[2]] <- attr(Y,'levels')
    return(PV)
  } else {
    if(distance=="ddeuclidean") {
      # Use data driven Euclidean distance
      for(i in 1:nr) {		
        # Add new observation NewX[i,] to Xtmp
        Xtmp <- rbind(X,NewX[i,])
        
        # Compute data driven euclidean distances
        for(m in 1:dimension) {
          Xtmp[,m] <- Xtmp[,m]/sqrt(var(Xtmp[,m]))
        }
        
        di <- as.matrix(dist(Xtmp))
        
        for(th in 1:L) {
          # Add NewX[i] temporarily to group th:
          thIndex <- c(which(Y==th),n+1)
          HH <- rep(0,n+1)
          HH[thIndex] <- 1
          nvec_tmp <- nvec
          nvec_tmp[th] <- nvec_tmp[th] + 1
          
          # Assign weights
          WXtmp <- matrix(0,n+1,n+1)
          for(j in thIndex) {
            WXtmp[j,order(c(di[j,c(1:(n+1))]))] <- W
          }
          
          # Compute PV[i,th]
          AA <- (WXtmp %*% HH)
          Tv <- AA[thIndex]
          PV[i,th] <- sum(Tv <= Tv[length(Tv)]) / nvec_tmp[th]
        }
      }
      dimnames(PV)[[2]] <- attr(Y,'levels')
      return(PV)
      
    } else {
    # Use Euclidean distance
      for(i in 1:nr) {		
        # Add new observation NewX[i,] to Xtmp
        Xtmp <- rbind(X,NewX[i,])
        di <- as.matrix(dist(Xtmp))
        for(th in 1:L) {
           # Add NewX[i] temporarily to group th:
          thIndex <- c(which(Y==th),n+1)
          HH <- rep(0,n+1)
          HH[thIndex] <- 1
          nvec_tmp <- nvec
          nvec_tmp[th] <- nvec_tmp[th] + 1
          
          # Assign weights
          WXtmp <- matrix(0,n+1,n+1)
          for(j in thIndex) {
            WXtmp[j,order(c(di[j,c(1:(n+1))]))] <- W
          }
          
          # Compute PV[i,th]
          AA <- (WXtmp %*% HH)
          Tv <- AA[thIndex]
          PV[i,th] <- sum(Tv <= Tv[length(Tv)]) / nvec_tmp[th]
        }
      }
      dimnames(PV)[[2]] <- attr(Y,'levels')
      return(PV)
    }
  }
}

