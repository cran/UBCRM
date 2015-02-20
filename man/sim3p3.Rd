\name{sim3p3}
\alias{sim3p3}

\title{
Simulation of one dose-escalation study with the classical 3+3 design
}

\description{
Given a true rates vector of DLT probabilities, the function simulate a 3+3 dose-escalation design.
}

\usage{
sim3p3(truerate, seed = NULL)
}

\arguments{
  \item{truerate}{A nlevel-length vector of true rates for the DLTs.
}

  \item{seed}{If not empty, the seed to use for random generation.
}
}

\value{
    \item{data}{Study data.}

    \item{mtd}{If reached, an integer representing the MTD level.}

    \item{lastdose}{An integer representing the last experimented dose.}
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


\seealso{
\code{\link{troisPtrois}}, \code{\link{ssim3p3}}
}

\examples{
# A 3-dose study with 10%, 20% and 30% of true rates for toxicity
sim3p3(c(0.1,0.2,0.3))
}