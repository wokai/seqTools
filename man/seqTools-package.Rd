\name{seqTools-package}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{seqTools-package}
\alias{seqTools}
%
\docType{package}
%
\title{SeqTools: Bioconductor package for analysis of FASTQ and fasta files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{ Analyze read length, phred scores and alphabeth frequency
    and DNA k-mers on uncompressed and compressed files.}
%
\details{
\tabular{ll}{
Package: \tab seqTools\cr
Type: \tab Package\cr
Version: \tab 0.99.31\cr
Date: \tab 2013-10-14\cr
License: \tab GPL-2\cr
Depends: \tab methods\cr
}
}
%
\author{Wolfgang Kaisers
Maintainer: Wolfgang Kaisers <kaisers@med.uni-duesseldorf.de>
}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and
    the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) Count DNA k-mer
countDnaKmers("ATAAATA", 2)
# B) Quality check on FASTQ file
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
fq <- fastqq("test_l6.fq")
plotPhredQuant(fq, 1)
}
%
\keyword{FASTQ}
\keyword{fasta}
