\name{troisPtrois}
\alias{troisPtrois}

\title{
Dose escalation with the 3+3 design
}

\description{
The function gives the next level to include patients following a 3+3 design. Needs an updated input dataframe with the CreData() structure.
}

\usage{
troisPtrois(data = data, lastdose)
}

\arguments{
  \item{data}{
Study dataframe with CreData() structure.
}

  \item{lastdose}{
Integer representing the last experimented dose level.
}

}

\value{
    \item{nextdose}{An integer representing the next recommended dose to experiment.}

    \item{mtd}{If reached, an integer representing the MTD.}
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
\code{\link{sim3p3}}, \code{\link{ssim3p3}}
}

\examples{
# Study initialization
data<- CreData(5,c("5 mg/m2","7 mg/m2","10 mg/m2","15 mg/m2","20 mg/m2"))
data

# Three patients are treated at the dose 1, without any observed DLT:
data<- updata(data,lastdose=1,npt=3,ndlt=0)
data

# 3+3 design
troisPtrois(data,lastdose=1)
}

