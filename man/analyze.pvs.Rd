\encoding{utf8}
\name{analyze.pvs}
\alias{analyze.pvs}
\title{
Analyze P-Values
}
\description{
Graphical displays and quantitative analyses of a matrix of p-values.
}
\usage{ analyze.pvs(pv, Y = NULL, alpha = 0.05, roc = TRUE, pvplot = TRUE, cex = 1) }
\arguments{
  \item{pv}{ matrix with p-values, e.g. output of \code{\link{cvpvs}} or \code{\link{pvs}}. }
  \item{Y}{ optional. Vector indicating the classes which the observations belong to. }
  \item{alpha}{ test level, i.e. 1 - confidence level. }
  \item{roc}{ logical. If \code{TRUE} and \code{Y} is not \code{NULL}, ROC curves are plotted.   }
  \item{pvplot}{ logical. If \code{TRUE} or \code{Y} is \code{NULL}, the p-values are displayed graphically. }
  \item{cex}{ A numerical value giving the amount by which plotting text should be magnified relative to the default. }
}
\details{
Displays the p-values graphically, i.e. it plots for each p-value a rectangle. The area of this rectangle is proportional to the the p-value. The rectangle is drawn blue if the p-value is greater than \code{alpha} and red otherwise. \cr
If \code{Y} is not \code{NULL}, i.e. the class memberships of the observations are known (e.g. cross-validated p-values), then additionally it plots the empirical ROC curves and prints some empirical conditional inclusion probabilities \eqn{I(b,\theta)} and/or pattern probabilities \eqn{P(b,S)}. Precisely, \eqn{I(b,\theta)} is the proportion of training observations of class \eqn{b} whose p-value for class \eqn{\theta} is greater than \eqn{\alpha}, while \eqn{P(b,S)} is the proportion of training observations of class \eqn{b} such that the \eqn{(1 - \alpha)}-prediction region equals \eqn{S}.
}
\value{
\item{T}{ Table containing empirical conditional inclusion and/or pattern probabilities for each class \eqn{b}. In case of \eqn{L = 2} or \eqn{L=3} classes, all patterns \eqn{S} are considered. In case of \eqn{L > 3}, all inclusion probabilities and some special patters \eqn{S} are considered. }
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
\seealso{ \code{\link{cvpvs}, \link{pvs}} }
\examples{
X <- iris[c(1:49, 51:99, 101:149), 1:4]
Y <- iris[c(1:49, 51:99, 101:149), 5]
NewX <- iris[c(50, 100, 150), 1:4]

cv <- cvpvs(X,Y)
analyze.pvs(cv,Y)

pv <- pvs(NewX, X, Y, method = 'k', k = 10)
analyze.pvs(pv)
}
\keyword{ classif }
