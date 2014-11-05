\name{fastqKmerSubsetLocs}
%
\alias{fastqKmerSubsetLocs}
%
\title{
fastqKmerSubsetLocs function: Counts for a given DNA k-mer subset position wise
from FASTQ files.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
Reads (compressed) FASTQ files and counts for given DNA k-mer subset for each
    position in sequence. The k-mer subset is given by a vector of k-mer
    indices. k-mer indices can be obtained from DNA k-mers with the function
    \code{kMerIndex}.}
\usage{
fastqKmerSubsetLocs(filenames, k=4, kIndex)
}
\arguments{
    \item{filenames}{\code{character}. Vector of fastqKmerSubsetLocs file
        names. Files can be gz compressed.}
    \item{k}{\code{integer}. Length of counted DNA k-mers.}
    \item{kIndex}{\code{integer}. Numeric values which represent indices of
        DNA-k mers.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Maximal allowed value for k is 12.}
%
\value{\code{list}. The length of the list equals the number of given
    filenames. Contains for each given file a matrix. Each matrix has one row
    for each given \code{kIndex} and an additional row with counts for all
    other DNA k-mers (labeled \code{other}). The number of columns equals the
    maximal sequence length in the FASTQ file.}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and
    the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{
Wolfgang Kaisers
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
k <- 4
kMers <- c("AAAA", "AACC", "AAGG")
kIdx <- kMerIndex(kMers)
res <- fastqKmerSubsetLocs("test_l6.fq", k, kIdx)
}
\keyword{fastqKmerSubsetLocs}
\keyword{kmer}
