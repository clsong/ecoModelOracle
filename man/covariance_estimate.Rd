% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{covariance_estimate}
\alias{covariance_estimate}
\title{Calculate Covariance Estimates for Growth and Death}
\usage{
covariance_estimate(data, species_variable, growth_function, death_function)
}
\arguments{
\item{data}{A data frame containing the necessary variables.}

\item{species_variable}{The name of the species variable within the data frame (provided as a string).}

\item{growth_function}{The name of the growth function (provided as a string).}

\item{death_function}{The name of the death function (provided as a string).}
}
\value{
A named list containing the estimated covariances:
  * cov_growth: Covariance between the growth variable and the species variable.
  * cov_death: Covariance between the death variable and the species variable.
}
\description{
Estimates covariances between growth/death variables and a specified species variable.
}
\examples{
covariance_estimate(
  example_data,
  species_variable = "prey",
  growth_function = "prey",
  death_function = "prey*predator + prey^2"
)
}
