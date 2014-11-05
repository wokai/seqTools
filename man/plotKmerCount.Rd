\name{plotKmerCount}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{plotKmerCount}
\alias{plotKmerCount-methods}
\alias{plotKmerCount,Fastqq-method}
%
\title{plotKmerCount:
    Creation of plots DNA for k-mer counts from Fastqq objects.}
%
\description{ The function creates plots from counted DNA k-mers
    from Fastqq objects.}
%
\usage{plotKmerCount(object,index,mxey,main="K-mer count",...)}
%
\arguments{
    \item{object}{\code{Fastqq}: Object which contains collected values from
        \code{nFiles} FASTQ files.}
    %
    \item{index}{\code{integer}(optional): Index of FASTQ file(s) for which
        data is plotted. When value is missing, k-mer counts for all contained
        data is plotted.}
    %
    \item{mxey}{\code{integer}(optional): Maximal value for y axis, given by
        power of 2 (when mxey=4, then maximal ylim value is 2^4 = 16). Allows
        overriding of automatic calculated values.}
    %
    \item{main}{\code{character}(optional): Caption text which printed into
        the output.}
    %
    \item{...}{Additional parameters which are passed down to the plot
        function.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Values for i must be in \{1,...,nFiles\}. The function shrinks the
    k-mer count table down to size of 4096 (k = 6) when k > 6 in order to
    limit the complexity of the plot.}
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
%
\note{The static size of the retured k-mer array is 4^k.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata",package="seqTools")
setwd(basedir)
#
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4,
    probeLabel=c("g4", "g5"))
#
plotKmerCount(fq)
plotKmerCount(fq,1)
plotKmerCount(fq, 1:2)
#
}
\keyword{plotKmerCount}
\keyword{kmer}
\seealso{Fastqq-class}
