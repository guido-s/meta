\name{smoking}
\docType{data}
\alias{smoking}
\alias{lungcancer}
\title{Smoking example}

\description{Meta-analyses on the effect of smoking on mortality risk.

 Data have been reconstructed based on the famous Smoking and Health
 Report to the Surgeon General (Bayne-Jones S et al., 1964).  Data
 sets can be used to evaluate the risk of smoking on overall mortality
 and lung-cancer deaths, respectively. The person time is attributed
 such that the rate ratios are equal to the reported mortality ratios
 implicitely assuming that the data have arisen from a homogeneous age
 group; more detailed information by age is not available from the
 report. Note, the group of "non-smokers" actually consists of all
 participants except those who are smokers of cigarettes
 only. Information on real non-smokers is not available from the
 published Smoking and Health Report.
}

\usage{
data(smoking)

data(lungcancer)
}

\format{A data frame with the following columns:
  \describe{
  \item{study}{Study label}
  \item{participants}{Total number of participants}
  \item{d.smokers}{Number of deaths in smokers' group}
  \item{py.smokers}{Person years at risk in smokers' group}
  \item{d.nonsmokers}{Number of deaths in non-smokers' group}
  \item{py.nonsmokers}{Person years at risk in non-smokers' group}
  }
}

\source{
  Bayne-Jones S et al. (1964),
  Smoking and Health: Report of the Advisory Committee to the Surgeon
  General of the United States. U-23 Department of Health, Education,
  and Welfare. Public Health Service Publication No. 1103.
  \url{http://profiles.nlm.nih.gov/ps/retrieve/ResourceMetadata/NNBBMQ}
} 

\seealso{\code{\link{metainc}}}

\examples{
data(smoking)

m1 <- metainc(d.smokers, py.smokers,
              d.nonsmokers, py.nonsmokers,
              data=smoking, studlab=study)
print(m1, digits=2)

data(lungcancer)

m2 <- metainc(d.smokers, py.smokers,
              d.nonsmokers, py.nonsmokers,
              data=lungcancer, studlab=study)
print(m2, digits=2)
}

\keyword{datasets}
