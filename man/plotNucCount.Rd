\name{plotNucCount}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{plotNucCount}
\alias{plotNucCount-methods}
\alias{plotNucCount,Fastqq-method}
%
\title{plotNucCount: Plots nucleotide counts from Fastqq objects.}
%
\description{The function creates plots from nucleotide counts
            from Fastqq objects.}
%
\usage{plotNucCount(object, nucs=16, maxx,...)}
%
\arguments{
    \item{object}{\code{Fastqq}: Object which contains collected values from 
        \code{nFiles} FASTQ files.}
    %
    \item{nucs}{\code{integer}(optional): Index of nucleotides for which data
        is plotted. When value is missing, k-mer counts for all contained data
        is plotted.}
    %
    \item{maxx}{\code{integer}(optional): When given, nucleotide counts are 
        plotted for the first maxx nucleotide positions. This option is used
        for displaying detailed plots from the first read nucleotide positions
        (which are sometimes not equally distributed).}
    %
    \item{...}{(currently unused).}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Values for i must be in \{1,...,nFiles\}. The nucs index encodes for
    IUPAC characters as shown in the following table.
\tabular{rrrrrrrrrrr}{
1 \tab A \tab | \tab 6  \tab R \tab | \tab 11 \tab M \tab | \tab 16 \tab N\cr
2 \tab C \tab | \tab 7  \tab Y \tab | \tab 12 \tab B \tab | \tab 17 \tab .\cr
3 \tab G \tab | \tab 8  \tab S \tab | \tab 13 \tab D \tab | \tab 18 \tab -\cr
4 \tab T \tab | \tab 9  \tab W \tab | \tab 14 \tab H \tab | \tab 19 \tab =\cr
5 \tab U \tab | \tab 10 \tab K \tab | \tab 15 \tab V \tab | \tab 20 \tab ''\cr
}
When count values for 'A' are to be plotted, 'nucs' must be =1. When count
    values for 'GC' are to be plotted, 'nucs' must be c(2,3).
}
%
\value{None.}
%
\references{Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
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
#
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4,
    probeLabel=c("g4", "g5"))
#
plotNucCount(fq)
plotNucCount(fq, 1)
plotNucCount(fq, 1:2)
#
}
\keyword{plotNucCount}
\keyword{kmer}
\seealso{Fastqq-class}
