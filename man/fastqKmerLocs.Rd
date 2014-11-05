\name{fastqKmerLocs}
%
\alias{fastqKmerLocs}
%
\title{
fastqKmerLocs function: Counts DNA k-mers position wise from FASTQ files.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Reads (compressed) FASTQ files and counts for DNA k-mers
    for each position in sequence.}
%
\usage{fastqKmerLocs(filenames, k=4)}
%
\arguments{
    \item{filenames}{Vector of FASTQ file names. Files can be gz compressed.}
    %
    \item{k}{Length of counted DNA k-mers.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Maximal allowed value for k is 12.}
\value{\code{list}. The length of the list equals the number of given
    filenames. 
    Contains for each given file a matrix with 4^k rows and 
    (maxSeqLen - k + 1) columns (maxSeqLen= maximum read length).
    The matrix contains for each k-mer and 
    k-mer-start position the counted values.}
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
\note{The static size of the retured k-mer array is 4^k.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
res <- fastqKmerLocs("test_l10_ATCGN.fq", k=2)
res <- fastqKmerLocs("test_l10_atcg.fq", k=2)
res <- fastqKmerLocs("test_l10_ATCGN.fq", k=2)
res <- fastqKmerLocs("test_l6_multi_line.fq", k=2)
}
\keyword{fastqKmerLocs}
\keyword{kmer}
