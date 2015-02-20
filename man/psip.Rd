\name{psip}
\alias{psip}
\alias{psit}
\alias{psil}

\title{
Dose-Toxicity modelisation functions
}

\description{
psip corresponds to the power model psy(s,a)=s^exp(a).
psit corresponds to the  hyperbolic tangent model psy(s,a)= ((tanh(s)+1)/2)**a.
psil corresponds to the  logistic model psy(s,a) = exp(3+a*s)/(1+exp(3+a*s)).
}

\usage{
psip(sgl,a)
psit(sgl,a)
psil(sgl,a)
}

\arguments{
  \item{sgl}{Dose level singleton.
}

  \item{a}{Parameter.
}
}

\value{
Numeric value of the computed function.
}

\references{
O'Quigley J., Pepe M., Fisher L. (1990). Continual Reassessment Method: a practical design for Phase I clinical trials in cancer. Biometrics 46, 33-48.

O'Quigley J., Shen LZ. (1996). Continual Reassessment Method: a likelihood approach. Biometrics 52, 673-684.

Paoletti X., Kramar A. (2009). A comparison of model choices for the Continual Reassessment Method in phase I cancer trials. Statistics in Medecine 28, 3012-3028.

Chamorey Emmanuel. (2009). Methodologie des essais de phase precoce en cancerologie: evolution des schemas et apport de la pharmacologie. These.

Garret-Mayer Elizabeth. (2006). The Continual Reassessment Method for dose-finding studies: a tutorial. Clinical Trials: 57-71.
}

\author{
Benjamin Esterni, Baboukar Mane. Unite de Biostatistique et de Methodologie, Institut Paoli-Calmettes, Marseille, France.
}