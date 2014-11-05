\name{countGenomeKmers}
%
\alias{countGenomeKmers}
%
\title{countGenomeKmers: Counting K-mers in DNA sequences.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Counts K-mers of DNA sequences inside a vector of DNA sequences.
    The k-mers are searched in a set of search windows, 
    which are defined by \code{start} and \code{width} parameter.
    From each position of the search window, a DNA k-mer is 
    identified on the right hand side on the given DNA sequence.
    Each value in the \code{start} vector defindes the left border
    of a search window.
    The size of the search window is given by the appropriate value in the 
    \code{width} vector. 
    The function is intended to count DNA k-mers in selected regions 
    (e.g. exons) on DNA chromosomes while respecting strand orientation.}
%
\usage{countGenomeKmers(dna, seqid, start, width, strand, k)}
%
\arguments{
    \item{dna}{\code{character}. Vector of DNA sequences. \code{dna} must 
        not contain other characters than "ATCGN".
        Capitalization does not matter.
        When a 'N' character is found, the current DNA k-mer is skipped.}
    %
    \item{seqid}{\code{numeric}. Vector of (1-based) values describing the
    index of the analyzed sequences inside the given dna vector.}
    %
    \item{start}{\code{numeric}. Vector of (1-based) start positions for 
        reading windows.}
    %
    \item{width}{\code{numeric}. Vector of window width values.}
    %
    \item{strand}{\code{factor or numeric}. First factor level (or numeric: 1)
        value will be interpreted as (+)-strand. For any other values, the
        reversed complement sequence will be counted (in left direction from
        start value).}
    %
    \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs.
        Only a single value is allowed (length(n) = 1!)}  
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The function returns a matrix.
    Each colum contains the motif-count values for one frame.
    Each row represents one DNA motif.
    The DNA sequence of the DNA motif is given as row.name.}
%
\value{\code{matrix}.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
sq <- "TTTTTCCCCGGGGAAAA"
seqid <- as.integer(c(1, 1))
start <- as.integer(c(6, 14))
width <- as.integer(c(4, 4))
strand <- as.integer(c(1, 0))
k <- 2
countGenomeKmers(sq, seqid, start, width, strand, k)
}
\keyword{countGenomeKmers}
