\name{HealthRWM}
\alias{HealthRWM}
\docType{data}
\title{
Health data from Greene (2012) applications.
}
\description{
The dataset is used in Greene (2012) and is taken from Riphahn, Wambach,
Million (2003).
}
\usage{data("HealthRWM")}
\format{
  A data frame with 27326 observations on the following 25 variables.
  \describe{
    \item{\code{ID}}{Person-identification number}
    \item{\code{female}}{Female=1; male=0}
    \item{\code{year}}{Calendar year of the observation}
    \item{\code{age}}{Age in years}
    \item{\code{hsat}}{Health satisfaction, coded 0 (low) to 10 (high)}
    \item{\code{handdum}}{Handicapped=1; otherwise=0}
    \item{\code{handper}}{Degree of handicap in percent (0 to 100)}
    \item{\code{hhninc}}{Household nominal monthly net income in German
      marks/10,000}
    \item{\code{hhkids}}{Children under age 16 in the household=1;
      otherwise=0}
    \item{\code{educ}}{Years of schooling}
    \item{\code{married}}{Married=1; otherwise=0}
    \item{\code{haupts}}{Highest schooling degree is Hauptschul
      degree=1; otherwise=0}    
    \item{\code{reals}}{Highest schooling degree is Realschul degree=1;
      otherwise=0}
    \item{\code{fachhs}}{Highest schooling degree is Polytechnical
      degree=1; otherwise=0}
    \item{\code{abitur}}{Highest schooling degree is Abitur=1;
      otherwise=0}    
    \item{\code{univ}}{Highest schooling degree is university degree=1;
      otherwise=0}
    \item{\code{working}}{Employed=1; otherwise=0}
    \item{\code{bluec}}{Blue-collar employee=1; otherwise=0}
    \item{\code{whitec}}{White-collar employee=1; otherwise=0}
    \item{\code{self}}{Self-employed=1; otherwise=0}
    \item{\code{beamt}}{Civil servant=1; otherwise=0}
    \item{\code{docvis}}{Number of doctor visits in last three months,}
    \item{\code{hospvis}}{Number of hospital visits in last calendar year,}
    \item{\code{public}}{Insured in public health insurance=1; otherwise=0}
    \item{\code{addon}}{Insured by add-on insurance=1; otherwise=0}
  }
}

\source{
On Greene (2012) online resources, and on the Journal of Applied
Econometrics website
(http://qed.econ.queensu.ca/jae/2003-v18.4/riphahn-wambach-million/).
}
\references{
Riphahn, R.T. and Wambach, A. and Million, A. (2003), \emph{Incentive
  Effects in the Demand for Health Care: A Bivariate Panel Count Data
  Estimation}, Journal of Applied Econometrics, Vol. 18, No. 4,
387--405.

Green, W.H.. (2012). \emph{Econometric Analysis, 7th edition}, Prentice Hall.  
}
\examples{
###### Example 13.7 of Greene (2012)
####################################

## Selecting the same data point and scaling income
##########
data(HealthRWM)
dat88 <- subset(HealthRWM, year==1988 & hhninc>0)
dat88$hhninc <- dat88$hhninc/10000

### A guess start
thet0 <- c(b0=log(mean(dat88$hhninc)),b1=0,b2=0,b3=0)

## Table 13.2 First column
g <- hhninc~exp(b0+b1*age+b2*educ+b3*female)
res0 <- nls(g, dat88, start=thet0, control=list(maxiter=100))
summary(res0)$coef

## Table 13.2 Second column
## Trying very hard to reproduce the results, 
## Who is right?
h1 <- ~age+educ+female
model1 <- momentModel(g, h1, thet0, vcov="MDS", data=dat88)
res1 <- gmmFit(model1, control=list(reltol=1e-10, abstol=1e-10))
summary(res1)@coef

## Table 13.2 third column (close enough)
## Here a sandwich vcov is required because it is not
## efficient GMM
h2 <- ~age+educ+female+hsat+married
model2 <- momentModel(g, h2, thet0, vcov="MDS", data=dat88)
res2 <- gmmFit(model2, type="onestep")
summary(res2, sandwich=TRUE)@coef

## Table 13.2 fourth column (Can't get closer than that)
res3 <- gmmFit(model2)
summary(res3)@coef

# Lets see what happens if we start on Greene solution
             
update(res3, theta0=c(b0=-1.61192, b1=.00092, b2=.04647, b3=-.01517))

## No...
}
\keyword{datasets}
