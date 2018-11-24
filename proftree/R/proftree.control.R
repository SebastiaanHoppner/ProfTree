proftree.control <- function(minbucket = 7L,
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
                             ...)
{
  # minbucket
  minbucket <- as.integer(minbucket)
  if (minbucket < 1L) {
    warning("parameter \"minbucket\" must be positive, default used instead")
    minbucket <- 7L
  }
  
  # minsplit
  minsplit <- as.integer(minsplit)
  if (minsplit < 2 * minbucket) {
    warning("parameter \"minsplit\" must be at least twice as large as \"minbucket\", changed")
    minsplit <- 2 * minbucket
  }
  
  # maxdepth
  maxdepth <- as.integer(maxdepth)
  if (maxdepth > 15) {
    warning("computations may take extremely long for \"maxdepth\" > 15 (or even be infeasible)")
  } else if (maxdepth < 1) {
    warning("parameter \"maxdepth\" must be equal or larger than 1, default used")
    maxdepth <- 9
  }
  
  # lambda
  lambda <- as.numeric(lambda)
  if (lambda < 0) {
    warning("parameter \"lambda\" must be equal or larger than 0, default used")
    lambda <- 1
  }
  
  # niterations
  niterations <- as.integer(niterations)
  if (niterations < 100L) {
    warning("computations may be unreliable for \"niterations\" < 100")
  }
  
  # miniterations
  miniterations <- as.integer(miniterations)
  if (miniterations < 100L) {
    warning("computations may be unreliable for \"miniterations\" < 100")
  } else if (miniterations > niterations) {
    stop("\"miniterations\" must be equal or smaller than \"niterations\"")
  }
  
  # errortol
  errortol <- as.numeric(errortol)
  if (errortol < 0 | errortol > 1) {
    warning("\"errortol\" must be between 0 and 1, default used")
    errortol <- 0.0005
  } else if (errortol > 0.0005) {
    warning("computations may be unreliable for \"errortol\" > 0.0005")
  }
  
  # ntrees
  ntrees <- as.integer(ntrees)
  if (ntrees < 2L) {
    warning("parameter \"ntrees\" must be equal or larger than 2, default used")
    ntrees <- 100
  } else if (ntrees < 10L) {
    warning("computations may be unreliable for \"ntrees\" < 10")
  }
  
  # elitismrange
  elitismrange <- as.integer(elitismrange)
  if (elitismrange < 2) {
    warning("\"elitismrange\" must be equal or larger than 2, default used")
    elitismrange <- max(c(ceiling(0.05*ntrees),2))
  }
  
  # alpha
  alpha <- as.numeric(alpha)
  if (alpha <= 1) {
    warning("parameter \"alpha\" must be strictly greater than 1, default used")
    alpha <- 6
  }
  
  # beta
  beta <- as.numeric(beta)
  if (beta <= 1) {
    warning("parameter \"beta\" must be strictly greater than 1, default used")
    beta <- 14
  }
  
  # clv
  clv <- as.numeric(clv)
  if (clv <= 0) {
    warning("parameter \"clv\" must be strictly positive, default used")
    clv <- 200
  }
  
  # d
  d <- as.numeric(d)
  if (d <= 0 | d >= clv) {
    warning("parameter \"d\" must be strictly positive and strictly smaller than \"clv\", default used")
    d <- 10
  }
  
  # f
  f <- as.numeric(f)
  if (f <= 0) {
    warning("parameter \"f\" must be strictly positive, default used")
    f <- 1
  }
  
  # seed
  if (is.null(seed)) {
    seed <- as.integer(runif(1, max = 2^16))
  }
  seed <- as.integer(seed)
  if (seed < -1L) {
    warning("parameter \"seed\" must be non-negative, default used instead")
    seed <- -1L
  }
  
  # verbose
  verbose <- as.logical(verbose)
  
  op <- c(pmutatemajor = 0.2,
          pmutateminor = 0.2,
          pcrossover   = 0.2,
          psplit       = 0.2,
          pprune       = 0.2)
  operatorprob <- unlist(operatorprob)
  if (is.null(names(operatorprob)) & length(operatorprob) == 5L) {
    names(operatorprob) <- names(op)
  }
  if (is.null(names(operatorprob))) {
    warning("incorrect specification of \"operatorprob\", default used")
    operatorprob <- op
  }
  op[names(operatorprob)] <- operatorprob
  operatorprob <- op/sum(op)*100
  
  return(list(minbucket = minbucket,
              minsplit  = minsplit,
              maxdepth  = maxdepth,
              lambda    = lambda,
              niterations   = niterations,
              miniterations = miniterations,
              errortol      = errortol,
              ntrees        = ntrees,
              elitismrange  = elitismrange,
              alpha = alpha,
              beta  = beta,
              clv   = clv,
              d     = d,
              f     = f,
              operatorprob = operatorprob,
              seed    = seed,
              verbose = verbose))
}
