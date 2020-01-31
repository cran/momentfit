\name{Griliches}
\alias{Griliches}
\docType{data}
\title{Return to Education Data
}
\description{
  Labour data on 758 young workers between 16 and 30 years hold. Each
  observation provides information on one individual at two points in
  time: in 1980 (variable with 80) and in the year given be the
  YEAR (variable without 80). 
}
\usage{data("Griliches")}
\format{
A data.frame with 758 observations and 20 variables.
\describe{
  \item{RNS, RNS80}{Dummy for residency in the southern states}
  \item{MRT, MRT80}{Dummy for marital status (1 if married)}
  \item{SMSA, SMSA80}{Dummy for residency in metropolitan areas}
  \item{MED}{Mother's education in years}
  \item{IQ}{IQ score}
  \item{KWW}{"Knowledge of the World of Work" test score}
  \item{Year}{The year of the first observation}
  \item{AGE, AGE80}{Age in years}
  \item{S, S80}{Completed years of schooling}
  \item{EXPR, EXPR80}{Experience in years}
  \item{TENURE, TENURE80}{Tenure im years}
  \item{LW, LW80}{log wage}
}}

 \source{
Online complements of Fumio Hayashi (2000)
}
\references{
Griliches, Z. (1976). Wages of Very Young Men.
\emph{Journal of Political Economy}, \bold{84}, S69--S85.

Blackburn, M. and Neumark, D. (1992). Unobserved Ability, Efficiency
Wages, and Interindustry Wage Differentials.
\emph{Quarterly Journal of Economics}, \bold{107}, 1421--1436.

Hayashi, F. (2000). \emph{Econometrics}, New Jersey: Princeton
University Press.
}
\keyword{datasets}
