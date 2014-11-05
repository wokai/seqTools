\name{gcContentMatrix}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{gcContentMatrix}
\alias{gcContentMatrix-methods}
\alias{gcContentMatrix,Fastqq-method}
%
\title{
gcContentMatrix: Returns matrix with read counts for GC content.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Returns a matrix with read counts.}
%
\usage{
gcContentMatrix(object)
}
\arguments{
  \item{object}{\code{Fastqq}: Object from wich data is copied.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The matrix contains one column for each FASTQ file. Rows labeled 
    from 0 to 100 which represents percent (\%) GC content. The matrix contains
    numbers of reads with the respective proportion of GC (Row 2 contains
    number of reads with 2\% GC content).}
%
\value{\code{matrix}.}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and 
    the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata",package="seqTools")
setwd(basedir)
#
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4,
    probeLabel=c("g4","g5"))
fqm<-gcContentMatrix(fq)
}
\seealso{gcContent}
\keyword{gcContentMatrix}
\keyword{kmer}
\keyword{fastqq}
