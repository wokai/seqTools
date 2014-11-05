\name{plotPhredQuant}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{plotPhredQuant}
\alias{plotPhredQuant-methods}
\alias{plotPhredQuant,Fastqq-method}
%
\title{ plotPhredQuant: Plots the position wise 10\%, 25 \%, 50\%, 75\% 
    and 90 \% quantiles of phred values.}
%
\description{ The function creates plots which describes the position wise 
    distribution of phred quantiles in single FASTQ files.}
%
\usage{plotPhredQuant(object, i, main, ...)}
%
\arguments{
    \item{object}{\code{Fastqq}: Object which contains collected values from
        \code{nFiles} FASTQ files.}
    %
    \item{i}{\code{integer}(optional): Index FASTQ file for which phred
        quantiles are plotted.}
    %
    \item{main}{\code{integer}(optional): The main title displayed on top of
        the plot. When missing, a standard text is printed.}
    %
    \item{...}{Other arguments which are passed to the internally called
        plot function.}
}
%
\value{None.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata",package="seqTools")
#
setwd(basedir)
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4,
        probeLabel=c("g4", "g5"))
#
plotPhredQuant(fq, 2)
# Same plot without x-axis
plotPhredQuant(fq,2, xaxt="n")
}
\keyword{plotPhredQuant}
\keyword{kmer}
\seealso{Fastqq-class}
