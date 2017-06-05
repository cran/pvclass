\encoding{utf8}
\name{cvpvs.gaussian}
\alias{cvpvs.gaussian}
\title{ Cross-Validated P-Values (Gaussian) }
\description{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. The p-values are based on a plug-in statistic for the standard Gaussian model. The latter means that the conditional distribution of \eqn{X}, given \eqn{Y=y}, is Gaussian with mean depending on \eqn{y} and a global covariance matrix.
}
\usage{
cvpvs.gaussian(X, Y, cova = c('standard', 'M', 'sym'))
}
\arguments{
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{cova}{ estimator for the covariance matrix: \cr
               'standard': standard estimator, \cr
               'M': M-estimator, \cr
               'sym': symmetrized M-estimator. }
}
\details{
Computes cross-validated nonparametric p-values for the potential class memberships of the training data. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using a plug-in statistic for the standard Gaussian model with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations.
}
\value{ 
\code{PV} is a matrix containing the cross-validated p-values. Precisely, for each feature vector \code{X[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}.
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
\code{\link{cvpvs}, \link{cvpvs.knn}, \link{cvpvs.wnn}, \link{cvpvs.logreg}}
}
\examples{
X <- iris[, 1:4]
Y <- iris[, 5]

cvpvs.gaussian(X, Y, cova = 'standard')
}
\keyword{ classif }
