\name{FST.SummaryStat.test}
\alias{FST.SummaryStat.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using summary statistics to test the association between an quantitative/dichotomous outcome variable and a gene by a score type test allowing for multiple functional annotation scores.
}
\description{
This function tests a specific gene using summary statistics (score vector and its covariance matrix)
}
\usage{
FST.SummaryStat.test(score,Sigma,Z,weights,B=5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{score}{
The score vector of length p, where p is the total number of genetic variables.
}
  \item{Sigma}{
The p*p covariance matrix of the score vector
}
  \item{Z}{
Functional annotation scores, an p*q matrix where p is the total number of genetic variables and q is the number of functional annotation scores. Note that the first column in Z should be all 1 if the users want the original weights of SKAT/burden test to be included.
}
  \item{weights}{
A numeric vector of weights for genetic variants (The length should be same as the number of genetic variants in the set.). These weights are usually based on minor allele frequencies.
}
  \item{B}{
Number of Bootstrap replicates. The default is 5000.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{p.value}{
P-value of the set based generalized score type test.
}
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
## FST.SummaryStat.test tests a region.
# Input: score (a score vector), Sigma (the covariance matrix of the score vector)

library(FSTpackage)

data(FST.example)
score<-FST.example$score;Sigma<-FST.example$Sigma;Z<-FST.example$Z;weights<-FST.example$weights

# test with 5000 bootstrap replicates
result<-FST.SummaryStat.test(score,Sigma,Z,weights,B=5000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{test}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
