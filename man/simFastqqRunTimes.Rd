\name{simFastqqRunTimes}
%
\alias{simFastqqRunTimes}
%
\title{simFastqqRunTimes: For given values of k and nSeq the function creates
    FASTQ files with simulated data, collects k-mer data with the fastqq
    function and reports the run times for the data collection.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{ For each combination of the parameters k and nSeq, 
    the function writes one FASTQ file and collects the data.
    The FASTQ files are equally structured: 
    Each read contains 17 randomly selected DNA 6-mers.
    Therefore the read-length is always 102.}
%
\usage{simFastqqRunTimes(k, nSeq, filedir=".")}
%
\arguments{
    \item{k}{\code{numeric}. k-mer sizes which are passed to \code{fastqq}.
        Default value is 2:15.}
    %
    \item{nSeq}{\code{numeric}. Number of simulated reads in FASTQ-file.
        Default value is (100, 1000, ..., 10000000).}
    %
    \item{filedir}{\code{character}. The output can be placed in a separate
        directory. When not existant, the function tries to create 'filedir'. 
        The function throws an error when writing is not permitted in the given
        directory (Could not open file ...).}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The FASTQ files contain the parameter settings inside their filename. 
    The files are created with 'writeSimFastq'.}
%
\value{\code{data.frame}. 
    The data frame has four columns: id, k, nSeq and runtime.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
\dontrun{
res <- simFastqqRunTimes(k=2:9, nSeq=100000)
plot(runtime~k,res,type="b")
}
}
%
\keyword{simFastqqRunTimes}
