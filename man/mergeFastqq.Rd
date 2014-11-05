\name{mergeFastqq}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{mergeFastqq}
\alias{mergeFastqq-methods}
\alias{mergeFastqq,Fastqq-method}
%
\title{mergeFastqq: Merges two Fastqq object into one.}
%
\description{The contents of two given Fastqq objects are merged together 
                into one resulting Fastqq object.}
%
\usage{mergeFastqq(lhs,rhs)}
\arguments{
  \item{lhs}{\code{Fastqq}.}
  \item{rhs}{\code{Fastqq}.}
}
\details{The data on all FASTQ files in the two incoming objects is merged
    together.
    The object has the same internal structure as if the data from all FASTQ
    files had been collected by a separate call of \code{fastqq} on the merged 
    FASTQ file names of the arguments.
    Duplicated probeLabel's are separated by adding of consecutive numbers as
    suffix to all probeLabel's.
    When \code{lhs} and \code{rhs} contain kmer-counts for different \code{k} 
    (getK), the function uses the \code{meltDownK} mechanism in order to
    equalize the \code{k} values.
    Therefore it is possible to compare samples which were counted with
    different \code{k} (i.e. k-mer resolution).}
%
\note{Note that the \code{meltDownK} mechanism is assotiated with a change of
DNA k-mer count values. See 'meltDownK' help (note) for more information.}
%
\value{S4 Object of class 'Fastqq'.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir<-system.file("extdata",package="seqTools")
setwd(basedir)
#
lhs<-fastqq("g4_l101_n100.fq.gz",k=4,"g4")
rhs<-fastqq("g5_l101_n100.fq.gz",k=4,"g5")
fq<-mergeFastqq(lhs,rhs)
}
\keyword{mergeFastqq}
