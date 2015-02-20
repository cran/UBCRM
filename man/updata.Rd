\name{updata}
\alias{updata}

\title{
Update the CRM dataframe after new patients' collected data
}

\description{
This function uptdates the CRM dataframe (result of the CreData routine) with new treated patients or observed DLTs.
}

\usage{
updata(data = data, lastdose, npt, ndlt)
}

\arguments{
  \item{data}{
Dataframe to be updated.
}

  \item{lastdose}{
Integer representing the dose to be updated.
}

  \item{npt}{
Number of new treated patients.
}

  \item{ndlt}{
Number of DLTs among the npt patients.
}
}

\value{
Updated dataframe.
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
\code{\link{CreData}}
}

\examples{
# Study initialization
data<- CreData(5,c("5 mg/m2","7 mg/m2","10 mg/m2","15 mg/m2","20 mg/m2"))
data

# Three patients are treated at the dose 1, without any observed DLT:
data<- updata(data,lastdose=1,npt=3,ndlt=0)
data
}
