\encoding{utf8}
\name{pvclass-package}
\alias{pvclass-package}
\alias{pvclass}
\docType{package}
\title{ P-Values for Classification }
\description{
Computes nonparametric p-values for the potential class memberships of new observations as well as cross-validated p-values for the training data. The p-values are based on permutation tests applied to an estimated Bayesian likelihood ratio, using a plug-in statistic for the Gaussian model, 'k nearest neighbors' or 'weighted nearest neighbors'. \cr
Additionally, it provides graphical displays and quantitative analyses of the p-values.
}
\details{
\tabular{ll}{
Package: \tab pvclass\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-07-29\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
Use \code{\link{cvpvs}} to compute cross-validated p-values, \code{\link{pvs}} to classify new observations and \code{\link{analyze.pvs}} to analyze the p-values.
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

