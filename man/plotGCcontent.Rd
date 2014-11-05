\name{plotGCcontent}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{plotGCcontent}
\alias{plotGCcontent-methods}
\alias{plotGCcontent,Fastqq-method}
%
\title{plotGCcontent: Plots the proportions of relative GC content for all
    FASTQ files.}
%
\description{ The function creates plots on proportions of relative GC content.
    For each FASTQ file from wich data is contained, 
    one separate line is plotted. 
    A value of 0.1 at the proportion of 40 says that 0.1 \% of the reads have
    40 \% GC content.}
%
\usage{plotGCcontent(object, main, ...)}
%
\arguments{
    \item{object}{\code{Fastqq}: Object which contains collected values from 
        \code{nFiles} FASTQ files.}
    %
    \item{main}{\code{integer}(optional): The main title displayed on top of 
        the plot. When missing, a standard text is printed.}
    %
    \item{...}{Other arguments which are passed to the internally called 
        plot function.}
}
%
\details{The area under each plotted line adds up to 1.}
%
\value{None.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4, 
    probeLabel=c("g4", "g5"))
#
plotGCcontent(fq)
}
\keyword{plotGCcontent}
\keyword{kmer}
\seealso{Fastqq-class}
