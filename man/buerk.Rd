\encoding{utf8}
\name{buerk}
\alias{buerk}
\docType{data}
\title{
Medical Dataset
}
\description{
This data set collected by Dr. Bürk at the university hospital in Lübeck contains data of 21556 surgeries in a certain time period (end of the nineties). Besides the mortality and the morbidity it contains 21 variables describing the condition of the patient and the surgery.
}
\usage{data(buerk)}
\format{
  A data frame with 21556 observations on the following 23 variables.
  \describe{
    \item{\code{age}}{Age in years}
    \item{\code{sex}}{Sex (1 = female, 0 = male)}
    \item{\code{asa}}{ASA-Score (American Society of Anesthesiologists), describes the physical condition on an ordinal scale: \cr
1 = A normal healthy patient \cr
2 = A patient with mild systemic disease \cr
3 = A patient with severe systemic disease \cr
4 = A patient with severe systemic disease that is a constant threat to life \cr
5 = A moribund patient who is not expected to survive without the operation \cr
6 = A declared brain-dead patient whose organs are being removed for donor purposes}
    \item{\code{rf_cer}}{Risk factor: cerebral (1 = yes, 0 = no)}
    \item{\code{rf_car}}{Risk factor: cardiovascular (1 = yes, 0 = no)}
    \item{\code{rf_pul}}{Risk factor: pulmonary (1 = yes, 0 = no)}
    \item{\code{rf_ren}}{Risk factor: renal (1 = yes, 0 = no)}
    \item{\code{rf_hep}}{Risk factor: hepatic (1 = yes, 0 = no)}
    \item{\code{rf_imu}}{Risk factor: immunological (1 = yes, 0 = no)}
    \item{\code{rf_metab}}{Risk factor: metabolic (1 = yes, 0 = no)}
    \item{\code{rf_noc}}{Risk factor: uncooperative, unreliable (1 = yes, 0 = no)}
    \item{\code{e_malig}}{Etiology: malignant (1 = yes, 0 = no)}
    \item{\code{e_vascu}}{Etiology: vascular (1 = yes, 0 = no)}
    \item{\code{antibio}}{Antibiotics therapy (1 = yes, 0 = no)}
    \item{\code{op}}{Surgery indicated (1 = yes, 0 = no)}
    \item{\code{opacute}}{Emergency operation (1 = yes, 0 = no)}
    \item{\code{optime}}{Surgery time in minutes}
    \item{\code{opsepsis}}{Septic surgery (1 = yes, 0 = no)}
    \item{\code{opskill}}{Expirienced surgeond, i.e. senior physician (1 = yes, 0 = no)}
    \item{\code{blood}}{Blood transfusion necessary (1 = yes, 0 = no)}
    \item{\code{icu}}{Intensive care necessary (1 = yes, 0 = no)}
    \item{\code{mortal}}{Mortality (1 = yes, 0 = no)}
    \item{\code{morb}}{Morbidity (1 = yes, 0 = no)}
  }
}
\source{
L. Dümbgen, B.-W. Igl and A. Munk (2008).
P-Values for Classification.
\emph{Electronic Journal of Statistics \bold{2}}, 468--493, available at \url{http://arxiv.org/abs/0801.2934}.
}
\keyword{datasets}
