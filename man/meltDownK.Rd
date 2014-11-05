\name{meltDownK}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{meltDownK}
\alias{meltDownK-methods}
\alias{meltDownK,Fastqq-method}
%
\title{
meltDownK: Condensing DNA k-mer count data to lower k-value
    (i.e. shorter DNA motifs).
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
Returns a copy of given object where DNA k-mer counts and first DNA k-mer 
    count table are reduced in size. }
%
\usage{
meltDownK(object, newK)
}
%
\arguments{
    \item{object}{\code{Fastqq}: Object from wich data is copied.}
    \item{newK}{\code{integer}: New value for k. Must be >=1 and <= old k.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The function sums all count values which belong the the new motif
    up. The new motif is the new-k sized prefix of the given k-mer motif.}
%
\note{The \code{meltDownK} mechanism is assotiated with a change of
DNA k-mer count values (by its accumulative character).
Also, count values from down-melted tables are not identical to directly
counted values for lower k. 
For example counting 'AAAA' with k=1 yields four 'A'. 
Counting 'AAAA' with k=2 yields three 'AA'.
As \code{meltDownK} sums up count values by prefix k-mers, the melted count
table for the second (k=2) count will return three 'A'.
Another source for differences may be N-nucleotides. Counting 'AANA' returns
three 'A' (using k=1) but only one 'AA' for k=2.}
%
\value{S4 Object of class 'Fastqq'.}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and
    the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
fq<-fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4, 
    probeLabel=c("g4","g5"))
fqm <- meltDownK(fq, 2)
}
\keyword{meltDownK}
\keyword{kmer}
\keyword{FASTQ}
\keyword{Fastqq}
