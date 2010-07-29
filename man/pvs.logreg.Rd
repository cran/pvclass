\encoding{utf8}
\name{pvs.logreg}
\alias{pvs.logreg}
\title{
P-Values to Classify New Observations (Penalized Multicategory Logistic Regression)
}
\description{
Computes nonparametric p-values for the potential class memberships of new observations. The p-values are based on 'penalized logistic regression'.
}
\usage{
pvs.logreg(NewX, X, Y, tau.o = 1,
           a0 = NULL, b0 = NULL,
           pen.method = c('vectors', 'simple', 'none'),
           progress = FALSE)
}
\arguments{
  \item{NewX}{ data matrix consisting of one or several new observations (row vectors) to be classified. }
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{tau.o}{ the penalty parameter (see section 'Details' below). }
  \item{a0, b0}{ optional starting values for logistic regression. }
  \item{pen.method}{ the method of penalization (see section 'Details' below). }
  \item{progress}{ optional parameter for reporting the status of the computations. }
}
\details{
Computes nonparametric p-values for the potential class memberships of new observations. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \code{Y[i]} equals \code{b}.
\cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using 'penalized logistic regression'. This means, the conditional probability of \eqn{Y = y}, given \eqn{X = x}, is assumed to be proportional to \eqn{exp(a_y + b_y^T x)}. The parameters \eqn{a_y}, \eqn{b_y} are estimated via penalized maximum log-likelihood. The penalization is either a weighted sum of the euclidean norms of the vectors \eqn{(b_1[j],b_2[j],\ldots,b_L[j])} (\code{pen.method=='vectors'}) or a weighted sum of all moduli \eqn{|b_{\theta}[j]|} (\code{pen.method=='simple'}). The weights are given by \code{tau.o} times the sample standard deviation (within groups) of the \eqn{j}-th components of the feature vectors. 
In case of \code{pen.method=='none'}, no penalization is used, but this option may be unstable.
}
\value{
\code{PV} is a matrix containing the p-values. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}.
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
\code{ \link{pvs}, \link{pvs.gaussian}, \link{pvs.knn}, \link{pvs.wnn}}
}
\examples{
X <- iris[c(1:49, 51:99, 101:149), 1:4]
Y <- iris[c(1:49, 51:99, 101:149), 5]
NewX <- iris[c(50, 100, 150), 1:4]

pvs.logreg(NewX, X, Y, tau.o=1, pen.method="vectors", progress=TRUE)

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

tmpi0 <- sample(1:n0.raw,size=3*n1,replace=FALSE)
tmpi1 <- sample(1:n1    ,size=  n1,replace=FALSE)

Xtrain <- rbind(X0[tmpi0[1:(n0-100)],],X1[1:(n1-100),])
Ytrain <- c(rep(1,n0-100),rep(2,n1-100))
Xtest <- rbind(X0[tmpi0[(n0-99):n0],],X1[(n1-99):n1,])
Ytest <- c(rep(1,100),rep(2,100))

PV <- pvs.logreg(Xtest,Xtrain,Ytrain,tau.o=2,progress=TRUE)
analyze.pvs(Y=Ytest,pv=PV,pvplot=FALSE)
}
}
\keyword{ classif }