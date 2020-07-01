\name{proftree.control}
\alias{proftree.control}
\title{Control for proftree}
\description{
  Various parameters that control aspects of the \code{proftree} fit.
}
\usage{
proftree.control(minbucket = 7L,
                 minsplit  = 20L,
                 maxdepth  = 9L,
                 lambda    = 0,
                 niterations   = 10000L,
                 miniterations = 1000L,
                 errortol      = 0.0005,
                 ntrees        = 100L,
                 elitismrange  = max(c(ceiling(0.05*ntrees),2)),
                 alpha = 6,
                 beta  = 14,
                 clv   = 200,
                 d     = 10,
                 f     = 1,
                 operatorprob = list(pmutatemajor = 0.2,
                                     pmutateminor = 0.2,
                                     pcrossover   = 0.2,
                                     psplit       = 0.2,
                                     pprune       = 0.2),
                 seed    = NULL,
                 verbose = FALSE,
                 \ldots)
}

\arguments{
  \item{minbucket}{the minimum sum of weights in a terminal node (default is 7). }
  \item{minsplit}{the minimum sum of weights in a node in order to be considered for splittin (default is 20). }
  \item{maxdepth}{maximum depth of the tree (default is 9). Note, that the memory
    requirements increase by the square of the m'aximum tree depth.}
  \item{lambda}{regulates the complexity part of the objective function.
    Increasing values of lambda encourage decreasing tree sizes (default is 0).}
  \item{niterations}{in case the run does not converge, it terminates
    after a specified number of iterations defined by \code{niterations} (default is 10000).}
  \item{miniterations}{minimum number of iterations before the algorithm terminates (default is 1000).}
  \item{errortol}{tolerance within the algorithm terminates (default is 0.0005).}
  \item{ntrees}{the number of trees in the population (default is 100).}
  \item{elitismrange}{integer; the algorithm terminates when the quality of the best \code{elitismrange} trees (default is 5\% of \code{ntrees}) stabilizes within a tolerance of \code{errortol} for 100 iterations, but not before \code{miniterations} iterations.}
  \item{alpha}{alpha parameter of unimodel beta distribution (default is 6).}
  \item{beta}{beta parameter of unimodel beta distribution (default is 14).}
  \item{clv}{constant CLV per retained customer (default is 200).}
  \item{d}{constant cost of retention offer (default is 10).}
  \item{f}{constant cost of contact (default is 1).}
  \item{operatorprob}{list or vector of probabilities for the selection
    of variation operators (default is 0.2 for each operator). May also be specified partially in which case
    the default values are still used for the unspecified arguments.
    Always scaled to sum to 100 percent.} 
  \item{seed}{a numeric seed to initialize the random number generator
    (for reproducibility). By default the seed is randomly drawn using
    \code{\link[stats]{runif}} in order to inherit the state of
    \code{\link[base]{.Random.seed}}. If set to \code{seed = -1L},    
    the random number generator is initialized by the system time.}
  \item{verbose}{should intermediate values of the objective function be shown (\code{TRUE}) or not (default is \code{FALSE})}
  \item{\dots}{additional arguments. }
}

\value{
  A list with the (potentially processed) control parameters.
}

