\name{updata}
\alias{updata}

\title{
update the crm dataframe after new patient's collected datas.
}

\description{
This function uptdates the crm dataframe (result of the CreData routine) with new treated patients or observed DLT's.
}

\usage{
updata(data = data, lastdose, npt, ndlt)
}

\arguments{
  \item{data}{
dataframe to be updated.
}

  \item{lastdose}{
integer representing the dose to be updated.
}

  \item{npt}{
number of new treated patients.
}

  \item{ndlt}{
number of DLT among the npt patients.
}
}

\value{
updated dataframe.
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
\code{\link{CreData}}
}

\examples{
# Study initialisation
data<- CreData(5,c("5 mg/m2","7 mg/m2","10 mg/m2","15 mg/m2","20 mg/m2"))
data

# three patients are treated at the dose 1, without any observed dlt:
data<- updata(data,lastdose=1,npt=3,ndlt=0)
data
}
