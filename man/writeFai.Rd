\name{writeFai}
%
\alias{writeFai}
%
\title{writeFai: Create FASTA index file.}
%
\description{The function reads a FASTA file and produces a FASTA index file
as output.}
%
\usage{writeFai(infiles, outfiles)}
%
\arguments{
    \item{infiles}{\code{character}. Vector of FASTA file names for which
        FASTA index is to be written.}
    %
    \item{outfiles}{\code{character}. Vector file names for writing FASTA
        index to.}
}
%
\value{None.}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and
    the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
\dontrun{
infile <- system.file("extdata", "small.fa", package="seqTools")
writeFai(infile, "small.fa.fai")
}
}
\keyword{writeFai}
