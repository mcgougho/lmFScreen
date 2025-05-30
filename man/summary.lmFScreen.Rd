% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.lmFScreen}
\alias{summary.lmFScreen}
\title{Summarize an lmFScreen Model}
\usage{
\method{summary}{lmFScreen}(object, ...)
}
\arguments{
\item{object}{An object of class lmFScreen containing model results.}

\item{...}{Additional arguments (currently unused).}
}
\value{
Invisibly returns a list containing:
\itemize{
\item Selective: A data frame with selective estimates, confidence intervals, and p-values
\item Standard: A data frame with standard estimates, confidence intervals, and p-values
}
}
\description{
Displays a structured summary of an lmFScreen model, including both selective
and naive estimates, confidence intervals, and p-values.
}
\details{
The output includes:
\itemize{
\item A header indicating the number of predictors
\item A table of selective estimates with confidence intervals and p-values
\item A table of naive estimates with confidence intervals and p-values
\item Significance stars are added to p-values:
\itemize{
\item for p < 0.05
** for p < 0.01
*** for p < 0.001
}
}
}
\examples{
data(mtcars)
result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
summary(result)

}
