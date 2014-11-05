\name{sim_fq}
%
\alias{sim_fq}
%
\title{sim_fq: Performs an experimental series of separation capabilities of
    hierarchical clustering (HC) based on DNA k-mers in FASTQ files using
    simulated DNA content.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{ Writes compressed FASTQ files where sequence sections
    contain concatenated k-mers which are uniformly distributed in the range
    of k-mers for given k.
    The function first writes a batch of randomly FASTQ files containing 
    randomly simulated DNA sequence.
    In a second step the function repeatedly writes FASTQ files with random
    DNA sequence where a fraction of the reads is 'contaminated' with given
    DNA k-mers.
    In a third step, for each set of simulated and contaminated files,
    a hierarchical cluster (HC) tree based on DNA k-mers is calculated.
    For each set of files, the size of the smaller fraction in the first 
    half of the tree is counted (perc).
    The value can be used as measure for separation capability of the
    HC algorithm.}
%
\usage{sim_fq(nRep=2, nContamVec=c(100, 1000), grSize=20, nSeq=1e4,
    k=6, kIndex=1365, pos=20)}
%
\arguments{
    \item{nRep}{\code{numeric}. Number of replicates for each combination of 
        each nContamVec value}
    %
    \item{nContamVec}{\code{numeric}. Vector with nContam (absolute number of
        contaminated reads) values.}
    %
    \item{grSize}{\code{numeric}. Number FASTQ files in control and 
        contamination group.}
    %
    \item{nSeq}{\code{numeric}. Number of reads per FASTQ file.}
    %
    \item{k}{\code{numeric}. k value used in fastqq function.}
    %
    \item{kIndex}{\code{numeric}. k-mer index of inserted k-mer(s). The k-mer
        index can be retreaved for a given k-mer with 'kMerIndex'. Default
        value is 1365 (="CCCCCC").}
    %
    \item{pos}{\code{numeric}. Determines at which position in sequence the
        k-mer is inserted. 1-based (1=first position).}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The function is intended to be used as explorative tool (not for
    routine quality assessment). There are some files written and there will
    be a lot of output on the terminal. It is therefore recommended to switch
    to a separate working directory and to run this function on a separate
    terminal. The function is not exported.}
%
\value{\code{data.frame} containing results of the counted perc values for
    each repetition of the simulation.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
kMerIndex("CCCCCC")
\dontrun{res <- seqTools:::sim_fq(nRep=2, nContamVec=c(10, 100),
                    grSize=4, nSeq=1e2)}
}
\keyword{sim_fq}
