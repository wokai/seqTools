\name{fastqq}
%
\alias{fastqq}
%
\title{
fastqq function: Reading summarizing information from FASTQ files.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
Reads read numbers, read lengths, counts per position alphabet frequencies,
    phred scores and counts per file DNA k-mers from (possibly compressed)
    FASTQ files.
}
\usage{
fastqq(filenames, k=6, probeLabel)
}
\arguments{
    \item{filenames}{Vector of FASTQ file names. Files can be gz compressed.}
    \item{k}{Length of counted DNA k-mers.}
    \item{probeLabel}{\code{character}: Textual label for each probe. When
    probeLabel and filenames have different length, a warning is thrown and
    the given labels are discarded.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{Maximal allowed value for k is 12.}
%
\value{S4 Object of class 'Fastqq'.}
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
fq <- fastqq("test_l6.fq")
fq <- fastqq("test_l6_multi_line.fq")
fq <- fastqq("non_exist.fq")
fq <- fastqq("test_l10_ATCGN.fq")
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4,
                                                    probeLabel=c("g4", "g5"))
}
\keyword{fastqq}
\keyword{Fastqq}
\keyword{kmer}
\seealso{Fastqq-class}
