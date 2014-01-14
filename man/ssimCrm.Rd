\name{ssimCrm}
\alias{ssimCrm}

\title{
Simulation of n dose-escalation study with the Continual Reassment Method.
}

\description{
The ssimCrm function simulates n dose escalation study with the CRM and provides summarized results
}

\usage{
ssimCrm(prior, n, firstdose = NA, truerate = prior, cohortsize = 3, target = 1/3,
nptmax = 24, nmaxmtd = nptmax, nmaxdose = nptmax, sd = 1.34, approach = "bayes",
method = "fpost", model = "power", nextlevel = "ntarget", upskipping = F,
downskipping = F, r = 2, seed = NULL)
}

\arguments{
  \item{prior}{numeric vector of prior DLTs probabilities
}

  \item{n}{Number of studies to simulate
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

  \item{method}{Estimation method for the posterior probabilities. "fpost" (default) estimates the mean of the posterior distribution of the paramater alpha (hat_alpha=E[alpha]) and uses it in psy(hat_alpha,...). "ppostp" and "pposts" directly estimate the mean of the posterior DLT probability. "ppostp" uses prior as singletons wheras "pposts" calculates appropriates singletons (see ail, ait or aip functions).
}

  \item{nextlevel}{character option used for determining the next dose level: "ntarget" (default) if the next level is choosed as the closest level to the desired target (may be higher than target). "utarget" if the next level is the closest level with the restriction to be lower than the target value.
}

  \item{upskipping}{boolean option used for determining the next dose level: If TRUE no level skip in escalation will be allowed. If FALSE (default) the level skips will be permitted.
}

  \item{downskipping}{boolean option used for determining the next dose level: If TRUE no level skip in desescalation will be allowed. If FALSE (default) the level skips will be permitted.
}

  \item{r}{Integer, number of digits for percentages in output
}

  \item{seed}{If non empty, the seed to use for random generation.
}
}

\value{
    \item{data}{summarized result in a "np1" view}

    \item{norecommendation}{Percentage of studies with no recommendation for the MTD (in case of the first level is considered as toxic}

    \item{mean.npt}{Mean number of enrolled patients}
    
    \item{mean.ndlt}{Mean number of observed DLTs}

    \item{mean.lastdose}{Mean last experimented dose level}
    
    \item{mean.prob}{Mean of posterior dlt probabilities}
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
\code{\link{Crm}}, \code{\link{simCrm}}
}

\examples{
ssimCrm(c(0.1,0.2,0.3,0.35,0.45),firstdose=1,target=0.33,n=100)
}