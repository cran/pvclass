analyze.pvs <-
function(pv, Y = NULL, alpha = 0.05,roc=TRUE,pvplot=TRUE){
  if(!is.null(Y)){
    Y <- factor(Y)
    T <- pvclass:::prtable(Y=unclass(Y),pv=pv,alpha=alpha)
    if(roc==TRUE & pvplot==TRUE){ par(ask=TRUE) }
    if(roc==TRUE){ pvclass:::rocs(Y=unclass(Y),pv=pv) }
    if(pvplot==TRUE){ pvclass:::pvplot(Y=Y,pv=pv,alpha=alpha) }
    invisible(T)
  } else{
    pvclass:::pvplot(pv=pv,alpha=alpha)
  }
}
  
