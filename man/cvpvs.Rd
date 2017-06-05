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
\code{\link{cvpvs.gaussian}, \link{cvpvs.knn}, \link{cvpvs.wnn}, \link{cvpvs.logreg}, \link{pvs}, \link{analyze.pvs}}
}
\examples{
X <- iris[,1:4]
Y <- iris[,5]

cvpvs(X,Y,method='k',k=10,distance='d')
}
\keyword{ classif }