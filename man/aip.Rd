\name{aip}
\alias{aip}
\alias{ail1}
\alias{ait1}
\alias{ail2}
\alias{ait2}

\title{
Functions to calculate the appropriate dose level singletons
}

\description{
Pool of functions to calculate dose level singletons values:
aip, ail2 and ait2 calculate sgl in order that E[psy] = prior,
ail1 and ait1 calculate sgl in order that psy(sgl,1) = prior.
}

\usage{
aip(p_prior, sd = 1.34)
ait1(p_prior, a=1)
ail1(p_prior, a=1)
ait2(p_prior)
ail2(p_prior)
}

\arguments{
  \item{p_prior}{Prior toxicity probability.
}

  \item{sd}{Standard deviation in case of normal distribution for the parameter.
}

  \item{a}{Rate in case of exponential distribution for the parameter.
}
}

\value{
Numeric length(p-prior)-vector.
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


