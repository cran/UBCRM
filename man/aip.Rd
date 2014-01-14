\name{aip}
\alias{aip}
\alias{ail1}
\alias{ait1}
\alias{ail2}
\alias{ait2}

\title{
Functions to calculate the appropriated dose level singletons.
}

\description{
Pool of functions to calculated dose level singletons values.
aip, ail2 and ait2 calculate sgl in order that E[psy] = prior
ail1 and ait1 calculate sgl in order that psy(sgl,1) = prior
}

\usage{
aip(p_prior, sd = 1.34)
ait1(p_prior, a=1)
ail1(p_prior, a=1)
ait2(p_prior)
ail2(p_prior)
}

\arguments{
  \item{p_prior}{prior toxicity probability
}

  \item{sd}{standard deviation in case of normal distribution for the parameter
}

  \item{a}{rate in case of exponential distribution for the parameter
}
}

\value{
numeric length(p-prior)-vector
}

\references{
O'Quigley, J., Pepe, M., Fisher, L (1990). Continual reassessment method: A practical design for Phase I clinical trials in cancer. Biometrics 46, 33-48.

O'Quigley J, Shen LZ. (1996). Continual Reassessment Method: a likelihood approach. Biometrics 52, 673-684

Paoletti, X., Kramar, A. (2009). A comparison of model choices for the continual reassessment method in phase I cancer trials. Statistics in Medecine 28, 3012-3028.

Chamorey Emmanuel. (2009). Methodologie des essais de phase precoce en cancerologie : evolution des schemas et apport de la pharmacologie. These

Garret-Mayer Elizabeth. (2006). The continual reassessment method for dose-finding studies : a tutorial. Clinical Trials: 57-71.
}

\author{
Benjamin Esterni, Baboukar Mane. Unite de Biostatistique et de Methodologie, Institut Paoli-Calmettes, Marseille, France.
}


