\name{UBCRM-package}
\alias{UBCRM-package}
\alias{UBCRM}
\docType{package}
\title{
UBCRM is a package containing functions to simulate and conduct dose escalation phase I studies
}
\description{
Two designs are implemented in the package: the classical 3+3 and the Continual Reassessment Method.
Simulations tools are also available to estimate the operating characteristics of the methods with several user-dependent options.
}
\details{
\tabular{ll}{
Package: \tab UBCRM\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2015-02-12\cr
License: \tab public \cr
}

}
\author{
Benjamin Esterni, Baboukar Mane. Unite de Biostatistique et de Methodologie, Institut Paoli-Calmettes, Marseille, France.
}
\references{
O'Quigley J., Pepe M., Fisher L. (1990). Continual Reassessment Method: a practical design for Phase I clinical trials in cancer. Biometrics 46, 33-48.

O'Quigley J., Shen LZ. (1996). Continual Reassessment Method: a likelihood approach. Biometrics 52, 673-684.

Paoletti X., Kramar A. (2009). A comparison of model choices for the Continual Reassessment Method in phase I cancer trials. Statistics in Medecine 28, 3012-3028.

Chamorey Emmanuel. (2009). Methodologie des essais de phase precoce en cancerologie: evolution des schemas et apport de la pharmacologie. These.

Garret-Mayer Elizabeth. (2006). The Continual Reassessment Method for dose-finding studies: a tutorial. Clinical Trials: 57-71.
}
\keyword{UBCRM}

\examples{

data<- CreData(4)
prior<-c(.05,.1,.15,.2)

# One study simulation
simCrm(prior=prior, firstdose = 2, truerate = prior, cohortsize = 3, target = 1/3,
nptmax = 18, nmaxmtd = 6, nmaxdose = 18, sd = 1.34, approach = "bayes", model = "power",
method = "fpost", nextlevel = "ntarget", upskipping = TRUE, downskipping = FALSE,
lastdose = NA, graphic = FALSE, seed = 20130110)

# N simulations with CRM
# Power model, no up skipping, start at dose 2
res1<- ssimCrm(prior=prior, 100, firstdose = 2, truerate = prior, cohortsize = 3,
target = 1/3, nptmax = 18, nmaxmtd = 6, nmaxdose = 18, sd = 1.34, approach = "bayes",
method = "fpost", model = "power", nextlevel = "ntarget", upskipping = TRUE,
downskipping = FALSE, r = 2, seed=20130110)
res1

# Simulations with 3+3 design
res2<- ssim3p3(truerate=prior, 100, r = 2, seed=20130110)
res2
}
