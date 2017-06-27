\name{FST.prelim}
\alias{FST.prelim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The preliminary data management for FST (functional score tests)
}
\description{
Before testing a specific gene using a score type test, this function does the preliminary data management, such as fitting the model under the null hypothesis.
}
\usage{
FST.prelim(Y, X=NULL, id=NULL, out_type="C")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The outcome variable, an n*1 matrix where n is the total number of observations
}
  \item{X}{
An n*d covariates matrix where d is the total number of covariates.
}
  \item{id}{
The subject id. This is used to match the genotype matrix. The default is NULL, where the a matched phenotype and genotype matrix is assumed.
}
  \item{out_type}{
Type of outcome variable. Can be either "C" for continuous or "D" for dichotomous. The default is "C".
}


}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  It returns a list used for function FST.test().
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(FSTpackage)

# Load data example
# Y: outcomes, n by 1 matrix where n is the total number of observations
# X: covariates, n by d matrix
# G: genotype matrix, n by p matrix where n is the total number of subjects
# Z: functional annotation matrix, p by q matrix

data(FST.example)
Y<-FST.example$Y;X<-FST.example$X;G<-FST.example$G

# Preliminary data management
result.prelim<-FST.prelim(Y,X=X)

}

\keyword{preliminary work}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
