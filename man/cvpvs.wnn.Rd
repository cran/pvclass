\encoding{utf8}
\name{cvpvs.wnn}
\alias{cvpvs.wnn}
\title{ Cross-Validated P-Values (Weighted Nearest Neighbors) }
\description{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. The p-values are based on 'weighted nearest-neighbors'.
}
\usage{
cvpvs.wnn(X, Y, wtype = c('linear', 'exponential'), W = NULL,
          tau = 0.3, distance = c('euclidean', 'ddeuclidean',
          'mahalanobis'), cova = c('standard', 'M', 'sym'))
}
\arguments{
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{wtype}{ type of the weight function (see section 'Details' below). }
  \item{W}{ vector of the (decreasing) weights (see section 'Details' below). }
  \item{tau}{ parameter of the weight function. If \code{tau} is a vector or \code{tau = NULL}, the program searches for the best \code{tau}. For more information see section 'Details'. }
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
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \code{Y[i]} equals \code{b}.
\cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using 'weighted nearest neighbors' with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations. \cr
The (decreasing) weights for the observations can be either indicated with a \eqn{n} dimensional vector \code{W} or (\code{if W = NULL}) one of the following weight functions can be used: \cr
linear: \deqn{W_i = \max(1-\frac{i}{n}/\tau,0),}{W[i] = max(1-(i/n)/\tau,0),}
exponential: \deqn{W_i = (1-\frac{i}{n})^\tau.}{W[i] = (1-i/n)^{\tau}.}
If \code{tau} is a vector, the program searches for the best \code{tau}. To determine the best \code{tau} for the p-value \code{PV[i,b]}, the class label of the training observation \eqn{X[i,]} is set temporarily to \code{b} and then for all training observations with \code{Y[j] != b} the sum of the weights of the observations belonging to class \code{b} is computed. Then the \code{tau} which minimizes the sum of these values is chosen. \cr
If \code{W = NULL} and \code{tau = NULL}, \code{tau} is set to \code{seq(0.1,0.9,0.1)} if \code{wtype = "l"} and to \code{c(1,5,10,20)} if \code{wtype = "e"}. 
}
\value{
\code{PV} is a matrix containing the cross-validated p-values. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
If \code{tau} is a vector or \code{NULL} (and \code{W = NULL}), \code{PV} has an attribute \code{"opt.tau"},  which is a matrix and \code{opt.tau[i,b]} is the best \code{tau} for observation \code{X[i,]} and class \code{b} (see section 'Details'). \code{"opt.tau"} is used to compute the p-values.
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
\code{\link{cvpvs}, \link{cvpvs.gaussian}, \link{cvpvs.knn}, \link{cvpvs.logreg}}
}
\examples{
X <- iris[, 1:4]
Y <- iris[, 5]

cvpvs.wnn(X, Y, wtype = 'l', tau = 0.5)
}
\keyword{ classif }
