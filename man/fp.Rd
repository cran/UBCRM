\name{fp}
\alias{fp}
\alias{ft}
\alias{fl}

\title{Density functions
}

\description{
Density functions for the model parameter.
fp(a,sd) is the normal density: 1/(sd*sqrt(2*pi))*exp(-(a^2)/(2*sd^2)).
ft and fl are the exponential density (with a fixed rate = 1): exp(-a).
}

\usage{
fp(a, sd)
ft(a)
fl(a)
}

\arguments{
  \item{a}{Parameter.
}
  \item{sd}{Standard deviation.
}
}

\value{
Numeric value of the computed density.
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


