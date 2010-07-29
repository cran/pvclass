\encoding{utf8}
\name{pvs.wnn}
\alias{pvs.wnn}
\title{
P-Values to Classify New Observations (Weighted Nearest Neighbors)
}
\description{
Computes nonparametric p-values for the potential class memberships of new observations. The p-values are based on 'weighted nearest-neighbors'.
}
\usage{
pvs.wnn(NewX, X, Y, wtype = c('linear', 'exponential'), W = NA,
        tau = 0.3, distance = c('euclidean', 'ddeuclidean',
        'mahalanobis'), cova = c('standard', 'M', 'sym'))
}
\arguments{
  \item{NewX}{ data matrix consisting of one or several new observations (row vectors) to be classified. }
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{wtype}{ type of the weight function (see section 'Details' below). }
  \item{W}{ vector of the (decreasing) weights (see section 'Details' below). }
  \item{tau}{ parameter of the weight function (see section 'Details' below). If \code{tau} is a vector, the program searches for the best \code{tau}. If \code{tau = NA}, \code{tau} is set to \code{seq(0.1,0.9,0.1)} if \code{wtype = "l"} and to \code{c(1,5,10,20,30,40,50)} if \code{wtype = "e"}. }
  \item{distance}{ the distance measure: \cr
                   'euclidean':    fixed Euclidean distance, \cr
                   'ddeuclidean':  data driven Euclidean distance (component-wise standardization), \cr
                   'mahalanobis':  Mahalanobis distance. }
  \item{cova}{ estimator for the covariance matrix: \cr
               'standard': standard estimator, \cr
               'M': M-estimator, \cr
               'sym': symmetrized M-estimator. }
}
\details{
Computes nonparametric p-values for the potential class memberships of new observations. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using 'weighted nearest neighbors' with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations. \cr
The (decreasing) weights for the observation can be either indicated with a \eqn{n} dimensional vector \code{W} or (\code{if W = NA}) one of the following weight functions can be used: \cr
linear: \deqn{W_i = \max((1-\frac{i}{n})/\tau,0)}{W[i] = max(1-(i/n)/\tau,0),}
exponential: \deqn{W_i = (1-\frac{i}{n})^\tau}{W[i] = (1-i/n)^{\tau}.}
For the exponential weight function \eqn{\tau} should be in \eqn{(0,1]} and for the linear weight function it should be greater than 1. \cr
If \code{tau} is a vector, the program searches for the best \code{tau}, i.e. for each new observation \code{NewX[i,]} and each class \code{b} the cross-validated p-values of the training data augmented with \code{(NewX[i,],b)} are computed and then the best \code{tau} is chosen in the sense that the sum of \code{PV[i,b]} for \code{Y[i] != b} is minimal. If \code{tau = NA}, \code{tau} is set to \code{seq(0.1,0.9,0.1)} if \code{wtype = "l"} and to \code{c(1,5,10,20,30,40,50)} if \code{wtype = "e"}. 
}
\value{
\code{PV} is a matrix containing the p-values. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
If \code{tau} is a vector or \code{NA}, \code{PV} has an attribute \code{"opt.tau"}, which is a matrix and \code{opt.tau[i,b]} is the best \code{tau} for observation \code{NewX[i,]} and class \code{b} in the sense that the sum of \code{PV[i,b]} for \code{Y[i] != b} is minimal. \code{opt.tau[i,b]} is used to compute the p-value for observation \code{NewX[i,]} and class \code{b}.
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
\code{ \link{pvs}, \link{pvs.gaussian}, \link{pvs.knn}, \link{pvs.logreg}}
}
\examples{
X <- iris[c(1:49, 51:99, 101:149), 1:4]
Y <- iris[c(1:49, 51:99, 101:149), 5]
NewX <- iris[c(50, 100, 150), 1:4]

pvs.wnn(NewX, X, Y, wtype = 'l', tau = 0.5)
}
\keyword{ classif }