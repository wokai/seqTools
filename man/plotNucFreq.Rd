\name{plotNucFreq}
%
\alias{plotNucFreq}
\alias{plotNucFreq-methods}
\alias{plotNucFreq,Fastqq-method}
%
\title{plotNucFreq: Plots the position wise relative nucleotide content for
    nucleotides 'A','C','G','T'.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The function creates plots on position wise relative 
    nucleotide content single FASTQ files.}
%
\usage{plotNucFreq(object, i, main, maxx, ...)}
%
\arguments{
    \item{object}{\code{Fastqq}: Object which contains collected values from 
        \code{nFiles} FASTQ files.}
    %
    \item{i}{\code{integer}(optional): Index FASTQ file for which nucleotide
        frequencies are plotted.}
    %
    \item{main}{\code{integer}(optional): The main title displayed on top of 
        the plot. When missing, a standard text is printed.}
    %
    \item{maxx}{\code{integer}(optional): Determines the maximum sequence 
        position for which counts are plotted. Small values (e.g. 15) allow 
        plotting the distribution on the first nucleotides at larger resolution
        (see reference).}
    %
    \item{...}{Other arguments which are passed to the internally called plot
        function.}
}
%
\value{None.}
%
\author{Wolfgang Kaisers}
%
\references{
Hansen KD, Brenner SE, Dudoit S.
Biases in Illumina transcriptome sequencing caused by random hexamer priming.
Nucleic Acids Research 2010 Vol.38 No.12 e131, doi: 10.1093/nar/gkq224
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
#
fq <- fastqq(c("g4_l101_n100.fq.gz", "g5_l101_n100.fq.gz"), k=4, 
    probeLabel=c("g4","g5"))
#
plotNucFreq(fq, 2)
# Same plot without x-axis
plotNucFreq(fq, 2, xaxt="n")
#
plotNucFreq(fq, 1, maxx=15)
}
\keyword{plotNucFreq}
\keyword{kmer}
\seealso{Fastqq-class}
