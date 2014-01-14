\name{CreData}
\alias{CreData}

\title{
Creates a CRM dataframe.
}
\description{
Creates a n-row summary dataframe indicating the number of treated patients and observed DLT at each of the n dose-levels. This is the dataframe's structure that will be needed in the different functions of the UBCRM package.
}
\usage{
CreData(ndose = 3, dosenames = paste("dose", 1:ndose, sep = " "))
}

\arguments{
  \item{ndose}{

number of dose levels.
}
  \item{dosenames}{
  a ndose-lenght character vector of labels for the dose levels.

}
}

\value{
A ndose * 3 dataframe containing
\item{dose}{integer value 1..ndose ordering the doses}
\item{npt}{integer count of the treated patients at dose i}
\item{ndlt}{integer count of the observed DLT at dose i}
}

\references{
O'Quigley, J., Pepe, M., Fisher, L (1990). Continual reassessment method: A practical design for Phase I clinical trials in cancer. Biometrics 46, 33-48.

O'Quigley J, Shen LZ. (1996). Continual Reassessment Method: a likelihood approach. Biometrics 52, 673-684

Paoletti, X., Kramar, A. (2009). A comparison of model choices for the continual reassessment method in phase I cancer trials. Statistics in Medecine 28, 3012-3028.

Chamorey Emmanuel. (2009). Methodologie des essais de phase precoce en cancerologie : evolution des schemas et apport de la pharmacologie. These

Garret-Mayer Elizabeth. (2006). The continual reassessment method for dose-finding studies : a tutorial. Clinical Trials: 57-71.
}

\author{
Baboukar Mane, Benjamin Esterni. Unite de Biostatistique et de Methodologie, Institut Paoli-Calmettes, Marseille, France.
}

\seealso{
\code{\link{updata}}
}

\examples{
data<- CreData(5,c("5 mg/m2","7 mg/m2","10 mg/m2","15 mg/m2","20 mg/m2"))
data
}

