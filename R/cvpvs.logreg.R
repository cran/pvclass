cvpvs.logreg <- function(X,Y,tau.o=1,
	pen.method=c("vectors","simple","none"),progress=TRUE)
{
	pen.method <- match.arg(pen.method)
	n <- dim(X)[1]
	d <- dim(X)[2]
	Y <- factor(Y)
	Y <- unclass(Y)
	L <- max(Y)
	X <- as.matrix(X)
	
	if (progress)
	{
		print('Computation of cross-validated p-values',
			quote=FALSE)
		print(paste('for ',as.character(n),
			' training observations.'),
			quote=FALSE)
		print('Preliminary log. regression:',
			quote=FALSE)
	}
	tmp <- pvclass:::penlogreg(X,Y,tau.o,
		pen.method=pen.method,progress=progress)
	a0 <- tmp$a
	b0 <- tmp$b
	
	PV <- matrix(1,nrow=n,ncol=L)
	for (i in 1:n)
	{
		if (progress)
		{
			print(paste('Observation no. ',as.character(i),' ...'),
				quote=FALSE)
		}
		NewX <- X[i,]
		Xr <- X[(1:n)!=i,]
		Yr <- Y[(1:n)!=i]
		PV[i,] <- pvclass:::pvs.logreg(NewX,Xr,Yr,tau.o=tau.o,
			pen.method=pen.method,a0=a0,b0=b0)
	}
	return(PV)
}