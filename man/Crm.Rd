\name{Crm}
\alias{Crm}

\title{
Dose-escalation with the Continual Reassessment Method
}

\description{
The function gives the next level to include patients following a model-based design. Needs an updated input dataframe with the CreData() structure.
}

\usage{
Crm(Dk, prior, target = 1/3, nptmax = 24, nmaxmtd = 6, nmaxdose = nptmax, sd = 1.34,
approach = "bayes", model = "power", method = "fpost", nextlevel = "ntarget",
upskipping = F, downskipping = F, lastdose = NA)
}

\arguments{
  \item{Dk}{
Study dataframe with CreData() structure.
}

  \item{prior}{Numeric vector of prior DLTs probabilities.
}

  \item{target}{Target used for the MTD determination.
}

  \item{nptmax}{Maximum number of patients to include in the study.
}

  \item{nmaxmtd}{Maximum number of patients to be treated at the designated MTD. Assign a high value (=nptmax) to avoid such a stopping rule.
}

  \item{nmaxdose}{Maximum number of patients to be treated at the same dose. Assign a high value (=nptmax) to avoid such a stopping rule.
}

  \item{sd}{Standard deviation used in case of a normal distribution assumption for the parameter.
}

  \item{approach}{Character indicating the estimation method: "bayes" (default value) for CRM or "mle" for CRML.
}

  \item{model}{Character indicating the dose-DLT relationship model: "power", "tangent" or "logistic". More informations in the details section.
}

  \item{method}{Estimation method for the posterior probabilities. "fpost" (default) estimates the mean of the posterior distribution of the parameter alpha (hat_alpha=E[alpha]) and uses it in psy(hat_alpha,...). "ppostp" and "pposts" directly estimate the mean of the posterior DLT probability. "ppostp" uses prior as singletons whereas "pposts" calculates appropriates singletons (see ail, ait or aip functions).
}

  \item{nextlevel}{Character option used for determining the next dose level. "ntarget" (default) if the next level is chosen as the closest level to the desired target (may be higher than target). "utarget" if the next level is the closest level with the restriction to be lower than the target value.
}

  \item{upskipping}{Boolean option used for determining the next dose level. If TRUE no level skip in escalation will be allowed. If FALSE (default) the level skips will be permitted.
}

  \item{downskipping}{Boolean option used for determining the next dose level. If TRUE no level skip in desescalation will be allowed. If FALSE (default) the level skips will be permitted.
}

  \item{lastdose}{Integer representing the last experimented dose level.
}
}

\details{
Details of the 3 dose-DLT relationship proposed models:
"power" for the power model psy(s,a)=s^exp(a),
"tangent" for the hyperbolic tangent model psy(s,a)= ((tanh(s)+1)/2)**a,
"logistic" for the logistic model psy(s,a) = exp(3+a*s)/(1+exp(3+a*s)).
Note: power and tangent models are equivalent after an appropriate transformation on the parameter.
}

\value{
\item{nextdose}{An integer representing the next recommended dose to experiment.}
\item{mtd}{If reached, an integer representing the MTD.}
\item{prob}{Posterior DLTs probabilities.}
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
\code{\link{simCrm}}, \code{\link{ssimCrm}}
}

\examples{
data<- CreData(5)
data<- updata(data,1,3,0)
data<- updata(data,2,3,1)
data<- updata(data,2,3,1)
data
Crm(data,prior=c(0.1,0.15,0.25,0.35,0.45),target=0.3,nextlevel="ntarget",nptmax=24,nmaxmtd=6)
data<- updata(data,3,3,2)
data
Crm(data,prior=c(0.1,0.15,0.25,0.35,0.45),target=0.3,nextlevel="ntarget",nptmax=24,nmaxmtd=6)
}
