\encoding{utf8}
\name{cvpvs.knn}
\alias{cvpvs.knn}
\title{ Cross-Validated P-Values (k Nearest Neighbors) }
\description{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. The p-values are based on 'k nearest neighbors'.
}
\usage{
cvpvs.knn(X, Y, k = NA, distance = c('euclidean', 'ddeuclidean',
          'mahalanobis'), cova = c('standard', 'M', 'sym'))
}
\arguments{
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{k}{ number of nearest neighbors. If \code{k} is a vector, the program searches for the best \code{k} in the sense that the sum of \code{PV[i,b]} for \code{Y[i] != b} is minimal. If \code{k = NA}, \code{k} is set to 1:ceiling(length(Y)/2). }
  \item{distance}{ the distance measure: \cr
                   "euclidean":    fixed Euclidean distance, \cr
                   "ddeuclidean":  data driven Euclidean distance (component-wise standardization), \cr
                   "mahalanobis":  Mahalanobis distance. }
  \item{cova}{ estimator for the covariance matrix: \cr
               'standard': standard estimator, \cr
               'M': M-estimator, \cr
               'sym': symmetrized M-estimator. }
}
\details{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using 'k nearest neighbors' with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations.
}
\value{ 
\code{PV} is a matrix containing the cross-validated p-values. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
If \code{k} is a vector or \code{NA}, \code{PV} has an attribute \code{"opt.k"}, which is the best \code{k} in the sense that the sum of \code{PV[i,b]} for \code{Y[i] != b} is minimal. \code{"opt.k"} is used to compute the p-values.
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
\code{\link{cvpvs}, \link{cvpvs.gaussian}, \link{cvpvs.wnn}, \link{cvpvs.logreg}}
}
\examples{
X <- iris[, 1:4]
Y <- iris[, 5]

cvpvs.knn(X, Y, k = c(5, 10, 15))
}
\keyword{ classif }
