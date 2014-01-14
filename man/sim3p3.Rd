\name{sim3p3}
\alias{sim3p3}

\title{
Simulation of one dose esclation study with the classical 3+3 design.
}

\description{
Given a true rates vector of dlt probabilities, the function simulate a 3+3 dose escalation design.
}

\usage{
sim3p3(truerate, seed = NULL)
}

\arguments{
  \item{truerate}{A nlevel-length vector of true rates for the DLTs.
}

  \item{seed}{If non empty, the seed to use for random generation.
}
}

\value{
    \item{data}{study data}

    \item{mtd}{if reached, an integer representing the MTD level}

    \item{lastdose}{an integer representing the last experimented dose}
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


\seealso{
\code{\link{troisPtrois}}, \code{\link{ssim3p3}}
}

\examples{
# A 3-doses study with 10%, 20% and 30% of true rates for toxicity
sim3p3(c(0.1,0.2,0.3))
}