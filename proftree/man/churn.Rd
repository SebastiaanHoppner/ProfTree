\name{churn}
\docType{data}
\alias{churn}
\title{Customer churn data}
\description{
Customer churn data set from a South Korean telecom operator
}
\usage{data(churn)}
\format{
A data frame containing 623 observations and 11 variables.
  \describe{
    \item{\code{region}}{Region where the customer lives.}
    \item{\code{prod_num}}{Number that identifies the customer’s product.}
    \item{\code{active_months}}{Time since the customer joined the operator (in months).}
    \item{\code{contract_period}}{Length of the contract period.}
    \item{\code{revenue_avg}}{Average revenue.}
    \item{\code{nonpay_period}}{How long did the customer not pay the bills?}
    \item{\code{overdue_amt}}{Amount that the customer is overdue.}
    \item{\code{count_disconnect}}{Number of times the service was disconnected.}
    \item{\code{count_complaint}}{Number of filed complaints.}
    \item{\code{autopay}}{Did the customer use the automatic payment option?}
    \item{\code{churn}}{Did the customer leave the telecom operator? (target response)}
  }


}
\keyword{datasets}
\examples{
data(churn)
str(churn)
head(churn)
summary(churn)
}
