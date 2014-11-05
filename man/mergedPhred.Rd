\name{mergedPhred}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{mergedPhred}
\alias{mergedPhred-methods}
\alias{mergedPhred,Fastqq-method}
%
\alias{mergedPhredQuantiles}
\alias{mergedPhredQuantiles-methods}
\alias{mergedPhredQuantiles,Fastqq-method}
%
\alias{plotMergedPhredQuant}
\alias{plotMergedPhredQuant-methods}
\alias{plotMergedPhredQuant,Fastqq-method}
%
\title{
mergedPhred functions: Retrieving and plotting of phred quantities from whole
    \code{Fastqq} objects.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
The \code{Fastqq} objects contain position-wise counted phred values. 
    The \code{mergedPhred} function adds the counted values for all FASTQ files
    together into a single matrix. The matrix then again contains position-wise 
    counted phred values. The \code{mergedPhredQuantiles} and 
    \code{plotMergedPhredQuant} are analogues to the \code{phredQuantiles} 
    and \code{plotPhredQuant} functions.}
%
\usage{
mergedPhred(object)
mergedPhredQuantiles(object, quantiles)
plotMergedPhredQuant(object, main, ...)
}
%
\arguments{
    \item{object}{\code{Fastqq}: Object which contains collected values from 
        \code{nFiles} FASTQ files.}
    \item{quantiles}{\code{numeric}: Vector of quantiles. All values must 
        be in [0,1].}
    \item{main}{\code{character}: String wich is used as figure caption.
        Passed internally to plot function.}
    \item{...}{Optional arguments which are passed to the \code{plot} function
        in \code{plotMergedPhredQuant}.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The function adds the phred values from all contained FASTQ data.}
%
\value{\code{mergedPhred} returns a \code{matrix} with 94 rows and 
    (maxSeqLen + 1) columns. \code{mergedPhredQuantiles} returns a
    \code{data.frame} with one row for each given quantile and max(seqLen(.))
    columns. \code{plotMergedPhredQuant} returns nothing.}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and
the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771\\
Ewing B, Green P
Base-calling of automated sequencer traces using phred. II. Error probabilities.
Genome Research 1998 Vol. 8 No. 3 186-194
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata",package="seqTools")
setwd(basedir)
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4, 
    probeLabel=c("g4", "g5"))
#
ph <- mergedPhred(fq)
ph[25:35, 1:15]
pq <- mergedPhredQuantiles(fq, c(0.25, 0.5, 0.75))
plotMergedPhredQuant(fq)
#
}
\keyword{mergedPhred}
\keyword{mergedPhredQuantiles}
\keyword{plotMergedPhredQuant}
