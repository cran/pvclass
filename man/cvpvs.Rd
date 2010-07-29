\encoding{utf8}
\name{cvpvs}
\alias{cvpvs}
\title{ Cross-Validated P-Values }
\description{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data.
}
\usage{
cvpvs(X, Y, method = c('gaussian','knn','wnn', 'logreg'), ...)
}
\arguments{
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{method}{
       one of the following methods: \cr
       'gaussian': plug-in statistic for the standard Gaussian model, \cr
       'knn': k nearest neighbors, \cr
       'wnn': weighted nearest neighbors, \cr
       'logreg': multicategory logistic regression with \eqn{l1}-penalization.}
  \item{\dots}{ further arguments depending on the method (see \code{\link{cvpvs.gaussian},} \cr
\code{\link{cvpvs.knn}, \link{cvpvs.wnn}, \link{cvpvs.logreg}}). }
}
\details{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using a plug-in statistic for the Gaussian model, 'k nearest neighbors', 'weighted nearest neighbors' or multicategory logistic regression with \eqn{l1}-penalization (see \code{\link{cvpvs.gaussian}, \link{cvpvs.knn}, \link{cvpvs.wnn}, \link{cvpvs.logreg}}) with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations.
}
\value{
\code{PV} is a matrix containing the cross-validated p-values. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}.
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
\code{\link{cvpvs.gaussian}, \link{cvpvs.knn}, \link{cvpvs.wnn}, \link{cvpvs.logreg}, \link{pvs}, \link{analyze.pvs}}
}
\examples{
X <- iris[,1:4]
Y <- iris[,5]

cvpvs(X,Y,method='k',k=10,distance='d')
}
\keyword{ classif }