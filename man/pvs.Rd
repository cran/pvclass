\encoding{utf8}
\name{pvs}
\alias{pvs}
\title{
P-Values to Classify New Observations 
}
\description{
Computes nonparametric p-values for the potential class memberships of new observations.
}
\usage{
pvs(NewX, X, Y, method = c('gaussian', 'knn', 'wnn', 'logreg'), ...)
}
\arguments{
  \item{NewX}{ data matrix consisting of one or several new observations (row vectors) to be classified. }
  \item{X}{ matrix containing training observations, where each observation is a row vector. }
  \item{Y}{ vector indicating the classes which the training observations belong to. }
  \item{method}{
       one of the following methods: \cr
       'gaussian': plug-in statistic for the standard Gaussian model, \cr
       'knn': k nearest neighbors, \cr
       'wnn': weighted nearest neighbors, \cr
       'logreg': multicategory logistic regression with \eqn{l1}-penalization.} 
  \item{\dots}{ further arguments depending on the method (see \code{\link{pvs.gaussian}, \link{pvs.knn}, \link{pvs.wnn}, \link{pvs.logreg}}). }
}
\details{
Computes nonparametric p-values for the potential class memberships of new observations. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}. \cr
This p-value is based on a permutation test applied to an estimated Bayesian likelihood ratio, using a plug-in statistic for the Gaussian model, 'k nearest neighbors', 'weighted nearest neighbors' or multicategory logistic regression with \eqn{l1}-penalization (see \code{\link{pvs.gaussian}, \link{pvs.knn}, \link{pvs.wnn}, \link{pvs.logreg}}) with estimated prior probabilities \eqn{N(b)/n}. Here \eqn{N(b)} is the number of observations of class \eqn{b} and \eqn{n} is the total number of observations.
}
\value{
\code{PV} is a matrix containing the p-values. Precisely, for each new observation \code{NewX[i,]} and each class \code{b} the number \code{PV[i,b]} is a p-value for the null hypothesis that \eqn{Y[i] = b}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{pvs.gaussian}, \link{pvs.knn}, \link{pvs.wnn}, \link{pvs.logreg}, \link{cvpvs}, \link{analyze.pvs}}
}
\examples{
X <- iris[c(1:49, 51:99, 101:149), 1:4]
Y <- iris[c(1:49, 51:99, 101:149), 5]
NewX <- iris[c(50, 100, 150), 1:4]

pvs(NewX, X, Y, method = 'k', k = 10)
}
\keyword{ classif }