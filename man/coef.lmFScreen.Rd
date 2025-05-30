% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{coef.lmFScreen}
\alias{coef.lmFScreen}
\title{Extract Coefficients from an lmFScreen Model}
\usage{
\method{coef}{lmFScreen}(object, ...)
}
\arguments{
\item{object}{An object of class lmFScreen.}

\item{...}{Currently unused.}
}
\value{
Invisibly returns a data frame containing:
\itemize{
\item Predictor: Names of the predictor variables
\item Estimate (Selective): Selective coefficient estimates
\item Estimate (Naive): Naive coefficient estimates
}
}
\description{
This S3 method prints and returns the selective and naive coefficient estimates
from an lmFScreen model.
}
\examples{
data(mtcars)
mod <- lmFScreen(mpg ~ wt + hp, data = mtcars)
coef(mod)


}
