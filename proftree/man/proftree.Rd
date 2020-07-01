\name{proftree}
\alias{proftree}

\title{Profit Driven Decision Trees for Churn Prediction}

\description{Method for fitting a profit driven decision tree for customer churn prediction by using an evolutionary algorithm to optimize a regularized Expected Maximum Profit measure for customer churn (EMPC).}

\usage{
proftree(formula, data, subset, na.action, weights, control = proftree.control(\dots), \dots)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit,
    no interactions should be used.}
  \item{data}{the data as a data frame.}
  \item{subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}. These input arguments are currently not available.}
  \item{weights}{optional integer vector of case weights. This input argument is currently not available.}
  \item{control}{a list of control arguments specified via
    \code{\link{proftree.control}}.}
  \item{\dots}{arguments passed to \code{\link{proftree.control}}.}
}

\details{
ProfTree (H0ppner et al., 2018) builds a profit driven decision tree for predicting customer churn. The Expected Maximum Profit measure for customer Churn (EMPC) has been proposed in order to select the most profitable churn model (Verbraken et al., 2013). ProfTree integrates the EMPC metric directly into the model construction and uses an evolutionary algorithm for learning profit driven decision trees.

The resulting trees can be summarized and visualized by the \code{\link[partykit]{print.constparty}},
and \code{\link[partykit]{plot.constparty}} methods provided by the \pkg{partykit} package.
Moreover, the \code{\link[partykit]{predict.party}} method can be used to compute fitted responses,
probabilities (for classification trees), and nodes.
}

\value{
  An object of class \code{\link[partykit]{party}}.
}

\references{
Hoppner, S., Stripling, E., Baesens, B., vanden Broucke, S., and Verdonck, T. (2018). Profit driven decision trees for churn prediction. \emph{European Journal of Operational Research. DOI: 10.1016/j.ejor.2018.11.072}
}

\examples{
library(proftree)
library(EMP)

data(churn)

Lambda <- 0.5

ProfTree <- proftree(formula = churn ~ .,
                     data = churn[1:100, ],
                     control = proftree.control(lambda = Lambda,
                                                seed = 2020,
                                                verbose = TRUE))
plot(ProfTree)
print(ProfTree$info)

scores <- predict(ProfTree, newdata = churn[1:100, ], type = "prob")[, 2]

EMP <- empChurn(scores = scores, classes = churn$churn[1:100])$EMP
print(EMP)
print(EMP - Lambda * width(ProfTree))
print(ProfTree$info$evalfun)
}

\keyword{tree}
