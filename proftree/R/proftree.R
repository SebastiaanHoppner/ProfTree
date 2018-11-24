proftree <- function(formula, data, subset, na.action, weights, control = proftree.control(...), ...) 
{
  timer_start <- proc.time()
  # original call
  ocall <- match.call()
  
  # set random seed
  if (!is.null(control$seed)) set.seed(control$seed)
  
  # build model.frame
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  # check terms
  mt <- attr(mf, "terms")
  if (any(attr(mt, "order") > 1L)) {
    warning("interaction terms are ignored")
  }
  
  # extract weights
  weights <- model.weights(mf)
  if (is.null(weights)) {
    weights <- rep.int(1L, nrow(mf))
  } else {
    if (!isTRUE(all.equal(as.vector(weights), round(as.vector(weights))))) {
      warning('weights must be integers')
    }
    weights <- as.integer(weights)
    if (!is.numeric(weights)) {
      stop('weights must be numeric')
    } else if (any(weights < 0)) {
      stop('all weights must be positive')
    } else if (any(weights == 0L)) {
      mf <- mf[weights > 0, , drop = FALSE]
      weights <- weights[weights > 0]
    }
  }
  mf[["(weights)"]] <- NULL
  
  # reorder and check all variables
  mf <- mf[, c(2:ncol(mf), 1L), drop = FALSE]
  if (length(unique(mf[,ncol(mf)])) == 1L) {
    stop('dependent variable has no variation')
  }
  mf <- mf[, sapply(mf, function(var) length(unique(var))) > 1L] # drop variables with only one level
  if (is.null(dim(mf))) {
    stop('independent variables are all constant')
  }
  
  nVariables <- ncol(mf)
  nInstances <- nrow(mf)
  prediction <- array(0L, nInstances)
  finalPerformance <- 0
  #fitnesscurve <- array(999999, control$niterations)
  
  control$method <- 1
  if (!is.factor(mf[,nVariables])) {
    stop('dependent variable is not a factor')
  } else if (length(levels(mf[,nVariables])) < 2L) {
    stop("dependent variable has only", length(levels(mf[,nVariables])), " level(s)")
  }
  
  for (i in 1:(nVariables-1)) {
    if (is.character(mf[,i])) {
      mf[,i] <- as.factor( mf[,i] )
      warning(paste("character variable", names(mf[i]), "was converted to a factor"))
    } else if (is.logical(mf[,i])) {
      mf[,i] <- as.factor(mf[,i])
    } else if (class(mf[,i])[1] == "Date") {
      stop(paste("variable type \"Date\" of variable \"", names(mf[i]), "\" is not supported" ), sep = '""')
    }
  }
  
  # calculate the maximum number of internal nodes a tree with size control$maxdepth can have
  maxNode <- 1 # number of internal nodes
  if (control$maxdepth > 1) {
    for(i in 1:(control$maxdepth-1) ) {
      maxNode <- maxNode*2+1
    }
  }
  
  if (2*control$minbucket > sum(weights)-1) {
    stop(paste("no split could be found \n \"minbucket\" must be smaller than half the weighted number of observations in the training data"))
  } else if (control$minsplit > sum(weights)-1) {
    stop(paste("no split could be found \n \"minsplit\" must be smaller than the weighted number of observations in the training data"))
  }
  
  # splitvariables and splitpoints
  splitV <- array(-999999, maxNode)
  splitP <- array(-999999, maxNode)
  
  # mf is transformed into a one dimensional vector such it can be passed to the C++ routines
  ndata <- array(-1, (nInstances*nVariables))
  # maxCat is the maximum number of categories of the variables
  maxCat <- 0
  # varType is negative for nominal variables and positive for ordered variables
  # the absolute value of varType[i] is the number of destinct values of variable i
  varType <- array(-1, nVariables)
  j <- 1
  k <- 1
  
  # calculation of varType and maxCat
  for (i in 1:nVariables ) {
    ndata[(1 + nInstances*(i-1)) : (nInstances[1] + nInstances*(i-1))] <- mf[,i]
    if (is.factor(mf[,i])) {
      if (is.ordered(mf[,i])) {
        varType[i] <- length(levels(mf[,i]))
        j <- j + j
      } else {
        varType[i] <- -length(levels(mf[,i]))
        k <- k + 1
        if (abs(varType[i]) > maxCat & i < nVariables) {
          maxCat <- abs(varType[i])
        }
      }
    } else {
      varType[i] <- length(unique(mf[,i]))
      j <- j +1
    }
  }
  # for categorical splits
  csplit <- array(2, maxCat*maxNode)
  
  #specification of function tree
  tree <- function(nInstances, nVariables, varType, ndata,
                   weights, prediction, finalPerformance, #fitnesscurve,
                   splitV, splitP, csplit, maxNode, control) {
    on.exit(.C("freememory", PACKAGE = "proftree"))
    out <- .C("tree", PACKAGE = "proftree",
              as.integer(nInstances),
              as.integer(nVariables),
              as.integer(varType),
              as.double(ndata),
              as.integer(weights),
              as.integer(prediction),
              as.double(finalPerformance),
              #as.double(fitnesscurve),
              as.integer(splitV),
              as.double(splitP),
              as.integer(csplit),
              as.integer(maxNode),
              as.integer(control$minbucket),
              as.integer(control$minsplit),
              as.integer(control$niterations),
              as.integer(control$miniterations),
              as.double(control$errortol),
              as.integer(control$ntrees),
              as.integer(control$elitismrange),
              as.integer(control$operatorprob[["pmutatemajor"]]),
              as.integer(control$operatorprob[["pmutateminor"]]),
              as.integer(control$operatorprob[["pcrossover"]]),
              as.integer(control$operatorprob[["psplit"]]),
              as.integer(control$operatorprob[["pprune"]]),
              as.integer(control$method),
              as.double(control$alpha),
              as.double(control$beta),
              as.double(control$clv),
              as.double(control$d),
              as.double(control$f),
              as.double(control$lambda),
              as.logical(control$verbose),
              as.double(control$seed)
    )
    return(out)
  }
  
  # Call of the tree function
  out <- tree(nInstances,
              nVariables,
              varType,
              ndata,
              weights,
              prediction,
              finalPerformance,
              #fitnesscurve,
              splitV,
              splitP,
              csplit,
              maxNode, 
              control)
  
  if (out[[8]][1] < 0L) { # CHANGE TO out[[9]]
    stop("no split could be found", call.=FALSE)
  }
  if (out[[14]] == control$niterations) { # CHANGE TO out[[15]]
    warning("algorithm did not converge, computations may be unreliable", call.=FALSE)
  }
  
  mtree <- list()
  mtree$varType    <- varType
  mtree$splitV     <- out[[8]]  # CHANGE TO out[[9]]
  mtree$splitP     <- out[[9]]  # CHANGE TO out[[10]]
  mtree$csplit     <- out[[10]] # CHANGE TO out[[11]]
  mtree$maxdepth   <- control$maxdepth
  mtree$weights    <- weights
  mtree$prediction <- out[[6]]+1
  mtree$maxCat     <- maxCat
  mtree$seed       <- control$seed
  init <- .initializeNode(mtree)
  node <- init[[1]]
  gid  <- init[[2]]
  
  prediction <- array(-999999, length(mtree$prediction))
  for(i in 1:length(gid)) {
    prediction[mtree$prediction == gid[i]] <- i
  }
  fitted <- data.frame(mtree$weights, prediction, mf[nVariables])
  names(fitted) <- c("(weights)", "(fitted)", "(response)")
  #fitnesscurve <- out[[8]][which(out[[8]] < 999999)]
  timer_end <- proc.time() - timer_start
  
  partyObject <- party(node = node,
                       data = mf,
                       fitted = fitted,
                       terms = mt,
                       info = list(method        = "proftree",
                                   version       = "0.0.2",
                                   niterations   = out[[14]], # CHANGE TO out[[15]]
                                   miniterations = out[[15]], # CHANGE TO out[[16]]
                                   errortol      = out[[16]], # CHANGE TO out[[17]]
                                   ntrees        = out[[17]], # CHANGE TO out[[18]]
                                   elitismrange  = out[[18]], # CHANGE TO out[[19]]
                                   lambda        = out[[30]], # CHANGE TO out[[31]]
                                   seed          = mtree$seed,
                                   #fitnesscurve  = fitnesscurve,
                                   evalfun       = out[[7]],
                                   runtime       = timer_end[3],
                                   call          = ocall))
  class(partyObject) <- c("constparty", "party")
  
  return(partyObject)
}
