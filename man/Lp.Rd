\name{Lp}
\alias{Lp}
\alias{Lt}
\alias{Ll}

\title{
Likelihood functions
}

\description{
Lp is the likelihood function for the power model psy(s,a)=s^exp(a).
Lt is the likelihood function for the hyperbolic tangent model psy(s,a)= ((tanh(s)+1)/2)**a.
Ll is the likelihood function for the logistic model psy(s,a) = exp(3+a*s)/(1+exp(3+a*s)).
}

\usage{
Lp(a, data, sgl)
Lt(a, data, sgl)
Ll(a, data, sgl)
}

\arguments{
  \item{a}{Parameter.
}

  \item{data}{CRM dataframe with a CreData() structure.
}

  \item{sgl}{Dose level singleton.
}
}

\details{
# Power model
Lp<- function(a, data, sgl){
npt<- data$npt
ndlt<- data$ndlt
sapply(a,FUN=function(a){prod((psip(sgl,a)^ndlt)*((1-psip(sgl,a))^(npt-ndlt)))})
}

# Hyperbolic tangent model
Lt<- function(a, data, sgl){
npt<- data$npt
ndlt<- data$ndlt
sapply(a,FUN=function(a){prod((psit(sgl,a)^ndlt)*((1-psit(sgl,a))^(npt-ndlt)))})
}

# Logistic model
Ll<- function(a, data, sgl){
npt<- data$npt
ndlt<- data$ndlt
sapply(a,FUN=function(a){prod((psil(sgl,a)^ndlt)*((1-psil(sgl,a))^(npt-ndlt)))})
}
}

\value{
Numeric value of the computed likelihood.
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


