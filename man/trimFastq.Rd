\name{trimFastq}
%
\alias{trimFastq}
%
\title{trimFastq: Performs sequence removal, trimming (fixed and quality based)
    and nucleotide masking on FASTQ files.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Fastq files sometimes need to be preprocessed before alignment.
    Three different mechanisms come into use here:
    Discarding whole reads, trimming sequences and masking nucleotides.
    This function performs all three mechanisms 
    together in one step.
    All reads with insufficient phred are discarded. 
    The reads can be trimmed ad each terminal side (on trim of fixed size and
    a trim based on quality thresholds).
}
%
\usage{trimFastq(infile, outfile="keep.fq.gz", discard="disc.fq.gz",
        qualDiscard=0, qualMask=0, fixTrimLeft=0,
        fixTrimRight=0, qualTrimLeft=0, qualTrimRight=0,
        qualMaskValue=78, minSeqLen=0)
}
%
\arguments{
    \item{infile}{\code{character}. Input FASTQ file. Only one infile is 
        allowed per function call.}
    %
    \item{outfile}{\code{character}. Output FASTQ file.}
    %
    \item{discard}{\code{character}. Output file in which discarded reads
        are written.}
    %
    \item{qualDiscard}{\code{numeric}. All reads which contain one or more
        phred scores <qualDiscard will be discarded (i.e. output to discard).}
    %
    \item{qualMask}{\code{numeric}. All nucleotides for which phred 
        score < qualMask will be overwritten with \code{qualMaskValue}.}
    %
    \item{fixTrimLeft}{\code{numeric}. Prefix of this size will be trimmed.}
    %
    \item{fixTrimRight}{\code{numeric}. Suffix of this size will be trimmed.}
    %
    \item{qualMaskValue}{\code{numeric}. ASCII replace value for masked
        nucleotides}
    %
    \item{qualTrimLeft}{\code{numeric}. Prefix where all phred scores 
        are < qualTrimLeft will be trimmed.}
    %
    \item{qualTrimRight}{\code{numeric}. Suffix where all phred scores
        are < qualTrimRight will be trimmed.}
    %
    \item{minSeqLen}{\code{numeric}. All reads where sequence length after 
        (fixed and quality based) trimming is <minSeqLen will be descarded 
        (i.e. output to discard).}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The function divides the input file into two outputs: 
    The output file (contains the accepted reads) and the discard file
    (contains the excluded reads).
    After trim operations, the function checks for remaining read 
    length. When the read length is smaller than minSeqLen, the read will be 
    discarded.}
%
\value{Numeric. A vector of length 2 which contains the number of reads which
    are written to output and to discard}
%
\references{
Ewing B, Green P
Base-calling of automated sequencer traces using phred. II. Error probabilities.
Genome Research 1998 Vol. 8 No. 3 186-194
}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
trimFastq("sim.fq.gz", qualDiscard=10, qualMask=15, fixTrimLeft=2,
    fixTrimRight=2, qualTrimLeft=28, qualTrimRight=30, minSeqLen=5)
}
\keyword{trimFastq}
