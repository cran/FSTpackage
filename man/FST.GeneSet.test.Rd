\name{FST.GeneSet.test}
\alias{FST.GeneSet.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test the association between an quantitative/dichotomous outcome variable and a large gene-set by a score type test allowing for multiple functional annotation scores.
}
\description{
Once the preliminary work is done using "FST.prelim()", this function tests a specifc gene.
}
\usage{
FST.GeneSet.test(result.prelim,G,Z,GeneSetID,Gsub.id=NULL,weights=NULL,
B=5000,impute.method='fixed')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result.prelim}{
The output of function "FST.prelim()"
}
  \item{G}{
Genetic variants in the target gene, an n*p matrix where n is the subject ID and p is the total number of genetic variables. Note that the number of rows in G should be same as the number of subjects. ***The column name should be the variable name, in order to be matched with the GeneSetID.
}
  \item{Z}{
Functional annotation scores, an p*q matrix where p is the total number of genetic variables and q is the number of functional annotation scores. Note that the first column in Z should be all 1 if the users want the original weights of SKAT/burden test to be included.
}
  \item{GeneSetID}{
A p*2 matrix indicating the genes in which the variables are located, where the first column is the genes' name and the second column is the variables' name.
}
  \item{Gsub.id}{
The subject id corresponding to the genotype matrix, an n dimensional vector. This is in order to match the phenotype and genotype matrix. The default is NULL, where the order is assumed to be matched with Y, X and time.
}
  \item{weights}{
A numeric vector of weights for genetic variants (The length should be same as the number of genetic variants in the set). These weights are usually based on minor allele frequencies. The default is NULL, where the beta(1,25) weights are applied.
}
  \item{B}{
Number of Bootstrap replicates. The default is 5000.
}
  \item{impute.method}{
Choose the imputation method when there is missing genotype. Can be "random", "fixed" or "bestguess". Given the estimated allele frequency, "random" simulates the genotype from binomial distribution; "fixed" uses the genotype expectation; "Best guess" uses the genotype with highest probability.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{n.marker}{
number of heterozygous SNPs in the SNP set.
}
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
## FST.prelim does the preliminary data management.
# Input: Y, X (covariates)
## FST.test tests a region.
# Input: G (genetic variants), Z (functional annotation scores) and result of FST.prelim

library(FSTpackage)

# Load data example
# Y: outcomes, n by 1 matrix where n is the total number of observations
# X: covariates, n by d matrix
# G: genotype matrix, n by p matrix where n is the total number of subjects
# Z: functional annotation matrix, p by q matrix

data(FST.example)
Y<-FST.example$Y;X<-FST.example$X;G<-FST.example$G;Z<-FST.example$Z;GeneSetID<-FST.example$GeneSetID

# Preliminary data management
result.prelim<-FST.prelim(Y,X=X,out_type='D')

# test with 5000 bootstrap replicates
result<-FST.GeneSet.test(result.prelim,G,Z,GeneSetID,B=5000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{test}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
