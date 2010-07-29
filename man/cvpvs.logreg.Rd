\encoding{utf8}
\name{cvpvs.logreg}
\alias{cvpvs.logreg}
\title{ Cross-Validated P-Values (Penalized Multicategory Logistic Regression) }
\description{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. The p-values are based on 'penalized logistic regression'.
}
\usage{
cvpvs.logreg(X, Y, tau.o = 1,
             pen.method = c('vectors', 'simple', 'none'),
             progress = TRUE)
}
\arguments{
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{tau.o}{ the penalty parameter (see section 'Details' below). }
  \item{pen.method}{ the method of penalization (see section 'Details' below). }
  \item{progress}{ optional parameter for reporting the status of the computations. }
}
\details{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \code{Y[i]} equals \code{b}, based on the remaining training observations.
\cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using 'penalized logistic regression'. This means, the conditional probability of \eqn{Y = y}, given \eqn{X = x}, is assumed to be proportional to \eqn{exp(a_y + b_y^T x)}. The parameters \eqn{a_y}, \eqn{b_y} are estimated via penalized maximum log-likelihood. The penalization is either a weighted sum of the euclidean norms of the vectors \eqn{(b_1[j],b_2[j],\ldots,b_L[j])} (\code{pen.method=='vectors'}) or a weighted sum of all moduli \eqn{|b_y[j]|} (\code{pen.method=='simple'}). The weights are given by \code{tau.o} times the sample standard deviation (within groups) of the \eqn{j}-th components of the feature vectors. 
In case of \code{pen.method=='none'}, no penalization is used, but this option may be unstable.
}
\value{
\code{PV} is a matrix containing the cross-validated p-values. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}, based on the remaining training observations.
}
\references{
L. Dümbgen, B.-W. Igl and A. Munk (2008).
P-Values for Classification.
\emph{Electronic Journal of Statistics \bold{2}}, 468--493, available at \url{http://arxiv.org/abs/0801.2934}.

N. Zumbrunnen (2009).
P-Values for Weighted Nearest-Neighbor Classifiers.
Master Thesis, IMSV, Univ. of Bern, available at \url{http://www.imsv.unibe.ch/unibe/philnat/imsv/content/e6079/e6261/e9754/zumbrunnen2009_eng.pdf}.
}
\author{
Niki Zumbrunnen \email{niki.zumbrunnen@stat.unibe.ch} \cr
Lutz Dümbgen \email{lutz.duembgen@stat.unibe.ch} \cr
\url{http://www.imsv.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}
}
\seealso{ 
\code{\link{cvpvs}, \link{cvpvs.gaussian}, \link{cvpvs.knn}, \link{cvpvs.wnn}}
}
\examples{

X <- iris[, 1:4]
Y <- iris[, 5]

cvpvs.logreg(X, Y, tau.o=1, pen.method="vectors",progress=TRUE)

# A bigger data example: Buerk's hospital data.
\dontrun{
data(buerk)
X.raw <- as.matrix(buerk[,1:21])
Y.raw <- buerk[,22]
n0.raw <- sum(1 - Y.raw)
n1 <- sum(Y.raw)
n0 <- 3*n1

X0 <- X.raw[Y.raw==0,]
X1 <- X.raw[Y.raw==1,]

tmpi0 <- sample(1:n0.raw,size=n0,replace=FALSE)
tmpi1 <- sample(1:n1    ,size=n1,replace=FALSE)

X <- rbind(X0[tmpi0,],X1)
Y <- c(rep(1,n0),rep(2,n1))

str(X)
str(Y)

PV <- cvpvs.logreg(X,Y,
	tau.o=5,pen.method="v",progress=TRUE)

analyze.pvs(Y=Y,pv=PV,pvplot=FALSE)
}

}
\keyword{ classif }
