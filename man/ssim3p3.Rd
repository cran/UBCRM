\name{ssim3p3}
\alias{ssim3p3}

\title{
Simulation of n dose-escalation study with the 3+3 design
}

\description{
The ssim3p3 function simulates n dose-escalation study with the 3+3 design and provides summarized results.
}

\usage{
ssim3p3(truerate, n, r = 2, seed = NULL)
}

\arguments{
  \item{truerate}{A nlevel-length vector of true rates for the DLTs.
}

  \item{n}{Number of studies to simulate.
}

  \item{r}{Integer, number of digits for percentages in output.
}

  \item{seed}{If not empty, the seed to use for random generation.
}
}

\value{
    \item{data}{Summarized result in a "np1" view.}

    \item{norecommendation}{Percentage of studies with no recommendation for the MTD (in case of the first level is considered as toxic).}

    \item{mean.npt}{Mean number of enrolled patients.}

    \item{mean.ndlt}{Mean number of observed DLTs.}

    \item{mean.lastdose}{Mean last experimented dose level.}
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
\code{\link{troisPtrois}}, \code{\link{sim3p3}}
}

\examples{
ssim3p3(c(0.1,0.2,0.25,0.35),100)
}
