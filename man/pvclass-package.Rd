\encoding{utf8}
\name{pvclass-package}
\alias{pvclass-package}
\alias{pvclass}
\docType{package}
\title{ P-Values for Classification }
\description{
Computes nonparametric p-values for the potential class memberships of new observations as well as cross-validated p-values for the training data. The p-values are based on permutation tests applied to an estimated Bayesian likelihood ratio, using a plug-in statistic for the Gaussian model, 'k nearest neighbors', 'weighted nearest neighbors' or 'penalized logistic regression'. \cr
Additionally, it provides graphical displays and quantitative analyses of the p-values.
}
\details{
Use \code{\link{cvpvs}} to compute cross-validated p-values, \code{\link{pvs}} to classify new observations and \code{\link{analyze.pvs}} to analyze the p-values.
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
\examples{
X <- iris[c(1:49, 51:99, 101:149), 1:4]
Y <- iris[c(1:49, 51:99, 101:149), 5]
NewX <- iris[c(50, 100, 150), 1:4]

cv <- cvpvs(X,Y)
analyze.pvs(cv,Y)

pv <- pvs(NewX, X, Y, method = 'k', k = 10)
analyze.pvs(pv)
}
\keyword{ package }
\keyword{ classif }

