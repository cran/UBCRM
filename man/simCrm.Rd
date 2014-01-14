\name{simCrm}
\alias{simCrm}

\title{
Simulation of one dose esclation study with the Continual Reassment Method.
}

\description{
Given prior and true rates vectors of dlt probabilities, the function simulate a CRM dose escalation design.
}

\usage{
simCrm(prior, firstdose = NA, truerate = prior, cohortsize = 3, target = 1/3,nptmax = 24,
nmaxmtd = nptmax, nmaxdose = nptmax, sd = 1.34, approach = "bayes", model = "power",
method = "fpost", nextlevel = "ntarget", upskipping = F, downskipping = F, lastdose = NA,
graphic = F, seed = NULL)
}

\arguments{
  \item{prior}{numeric vector of prior DLTs probabilities
}

  \item{firstdose}{Integer representing the dose at which the first cohort will be treated.
}

  \item{truerate}{A nlevel-length vector of true rates for the DLTs.
}

  \item{cohortsize}{Size of the cohort size. Default value = 3
}

  \item{target}{Target used for the MTD determination
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

  \item{model}{character indicating the dose-dlt relationship model: "power", "tangent", or "logistic".
}

  \item{method}{stimation method for the posterior probabilities. "fpost" (default) estimates the mean of the posterior distribution of the paramater alpha (hat_alpha=E[alpha]) and uses it in psy(hat_alpha,...). "ppostp" and "pposts" directly estimate the mean of the posterior DLT probability. "ppostp" uses prior as singletons wheras "pposts" calculates appropriates singletons (see ail, ait or aip functions).
}

  \item{nextlevel}{character option used for determining the next dose level: "ntarget" (default) if the next level is choosed as the closest level to the desired target (may be higher than target). "utarget" if the next level is the closest level with the restriction to be lower than the target value.
}

  \item{upskipping}{boolean option used for determining the next dose level: If TRUE no level skip in escalation will be allowed. If FALSE (default) the level skips will be permitted.
}

  \item{downskipping}{boolean option used for determining the next dose level: If TRUE no level skip in desescalation will be allowed. If FALSE (default) the level skips will be permitted.
}

  \item{lastdose}{Integer representing the last experimented dose level.
}

  \item{graphic}{bolean option for graphic generation.
}

  \item{seed}{If non empty, the seed to use for random generation.
}
}

\value{
    \item{data}{study data}

    \item{dose}{integer vector representing for each cohort the experimented dose levels}

    \item{ndlt}{integer vector representing for each cohort the number of observed DLTs}

    \item{mtd}{if reached, an integer representing the MTD level}

    \item{lastdose}{an integer representing the last experimented dose}

    \item{prob}{posterior dlt probabilities}
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
\code{\link{Crm}}, \code{\link{ssimCrm}}
}

\examples{
simCrm(c(0.1,0.2,0.3,0.35,0.45),firstdose=1,target=0.33)
}