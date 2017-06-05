\encoding{utf8}
\name{pvs.knn}
\alias{pvs.knn}
\title{
P-Values to Classify New Observations (k Nearest Neighbors)
}
\description{
Computes nonparametric p-values for the potential class memberships of new observations. The p-values are based on 'k nearest neighbors'.
}
\usage{
pvs.knn(NewX, X, Y, k = NULL, distance = c('euclidean', 'ddeuclidean',
        'mahalanobis'), cova = c('standard', 'M', 'sym'))
}
\arguments{
  \item{NewX}{ data matrix consisting of one or several new observations (row vectors) to be classified. }
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{k}{ number of nearest neighbors. If \code{k} is a vector or \code{k = NULL}, the program searches for the best \code{k}. For more information see section 'Details'. }
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
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using 'k nearest neighbors' with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations. \cr
If \code{k} is a vector, the program searches for the best \code{k}. To determine the best \code{k} for the p-value \code{PV[i,b]}, the new observation \code{NewX[i,]} is added to the training data with class label \code{b} and then for all training observations with \code{Y[j] != b} the proportion of the \code{k} nearest neighbors of \code{X[j,]} belonging to class \code{b} is computed. Then the \code{k} which minimizes the sum of these values is chosen. \cr
If \code{k = NULL}, it is set to 2:ceiling(length(Y)/2).
}
\value{
\code{PV} is a matrix containing the p-values. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
If \code{k} is a vector or \code{NULL}, \code{PV} has an attribute \code{"opt.k"}, which is a matrix and \code{opt.k[i,b]} is the best \code{k} for observation \code{NewX[i,]} and class \code{b} (see section 'Details'). \code{opt.k[i,b]} is used to compute the p-value for observation \code{NewX[i,]} and class \code{b}.
}
\references{
Zumbrunnen N. and Dümbgen L. (2017)
pvclass: An R Package for p Values for Classification.
\emph{Journal of Statistical Software \bold{78(4)}}, 1--19.
doi:10.18637/jss.v078.i04

Dümbgen L., Igl B.-W. and Munk A. (2008)
P-Values for Classification.
\emph{Electronic Journal of Statistics \bold{2}}, 468--493, available at \url{http://dx.doi.org/10.1214/08-EJS245}.

Zumbrunnen N. (2014)
P-Values for Classification – Computational Aspects and Asymptotics.
Ph.D. thesis, University of Bern, available at \url{http://boris.unibe.ch/id/eprint/53585}.
}
\author{
Niki Zumbrunnen \email{niki.zumbrunnen@gmail.com} \cr
Lutz Dümbgen \email{lutz.duembgen@stat.unibe.ch} \cr
\url{www.imsv.unibe.ch/duembgen/index_ger.html}
}
\seealso{ 
\code{\link{pvs}, \link{pvs.gaussian}, \link{pvs.wnn}, \link{pvs.logreg}}
}
\examples{
X <- iris[c(1:49, 51:99, 101:149), 1:4]
Y <- iris[c(1:49, 51:99, 101:149), 5]
NewX <- iris[c(50, 100, 150), 1:4]

pvs.knn(NewX, X, Y, k = c(5, 10, 15))
}
\keyword{ classif }