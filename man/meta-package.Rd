\name{meta-package}

\alias{meta-package}

\title{meta: Brief overview of methods and general hints}

\description{
R package \bold{meta} is a user-friendly general package providing
standard methods for meta-analysis and supporting Schwarzer et
al. (2015), \url{http://meta-analysis-with-r.org/}.
}

\details{
  R package \bold{meta} (Schwarzer, 2007) provides the following
  meta-analysis methods:
  \itemize{
  \item fixed effect and random effects meta-analysis (functions
    \code{\link{metabin}}, \code{\link{metacont}},
    \code{\link{metacor}}, \code{\link{metagen}},
    \code{\link{metainc}}, \code{\link{metaprop}}, and
    \code{\link{metarate}});
  \item several plots (\code{\link{forest}}, \code{\link{funnel}},
    Galbraith / \code{\link{radial}}, \code{\link{labbe}},
    \code{\link{baujat}}, \code{\link{bubble}});
  \item statistical tests (\code{\link{metabias}}) and trim-and-fill
    method (\code{\link{trimfill}}) to evaluate bias in meta-analysis;
  \item import data from 'RevMan 5' (\code{\link{read.rm5}}; see also
    \code{\link{metacr}});
  \item prediction interval, Hartung-Knapp and Paule-Mandel method for
    random effects model (arguments in meta-analysis functions);
  \item cumulative meta-analysis (\code{\link{metacum}}) and
    leave-one-out meta-analysis (\code{\link{metainf}});
  \item meta-regression (\code{\link{metareg}}; if R package
    \bold{metafor} is installed);
  \item generalised linear mixed models (\code{\link{metabin}},
    \code{\link{metainc}}, \code{\link{metaprop}}, and
    \code{\link{metarate}}; if R packages \bold{metafor}, \bold{lme4},
    \bold{numDeriv}, and \bold{BiasedUrn} are installed).
  }
  
  The following more advanced statistical methods are provided by
  add-on R packages:
  \itemize{
  \item frequentist methods for network meta-analysis (R package
    \bold{netmeta});
  \item advanced methods to model and adjust for bias in meta-analysis
    (R package \bold{metasens}).
  }
  
  See \code{\link{settings.meta}} to learn how to print and specify
  default meta-analysis methods used during your R session. For
  example, the following R command specifies the use of the
  Hartung-Knapp and Paule-Mandel methods, and the printing of
  prediction intervals for any meta-analysis in the current R session:
  \itemize{
  \item \code{settings.meta(hakn=TRUE, method.tau="PM", prediction=TRUE)}
  }

  Type \code{help(package = "meta")} for a listing of R functions
  available in \bold{meta}.

  Schwarzer (2007) is the preferred citation in publications for
  \bold{meta}. Type \code{citation("meta")} for a BibTeX entry of this
  publication.

  To report problems and bugs
  \itemize{
  \item type \code{bug.report(package = "meta")} if you do not use
    RStudio,
  \item send an email to Guido Schwarzer
    \email{sc@imbi.uni-freiburg.de} if you use RStudio.
  }

  The development version of \bold{meta} is available on GitHub
  \url{https://github.com/guido-s/meta}.
}

\references{
  Schwarzer G (2007),
  meta: An R package for meta-analysis.
  \emph{R News}, \bold{7}(3), 40--5.
  \url{https://cran.r-project.org/doc/Rnews/Rnews_2007-3.pdf}

  Schwarzer G, Carpenter JR and Rücker G (2015),
  \emph{Meta-Analysis with R (Use-R!)}.
  Springer International Publishing, Switzerland.
  \url{http://www.springer.com/gp/book/9783319214153}
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
