\name{countDnaKmers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{countDnaKmers}
%
\title{countDnaKmers: Counting k-mers in DNA sequence.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Counts occurrence of DNA k-mers in given DNA sequence. 
The k-mers are searched in a set of search windows,
which are defined by \code{start} and \code{width} parameter.
From each position of the search window, a DNA k-mer is identified
on the right hand side on the given DNA sequence. 
Each value in the \code{start} vector defines the left border
of a search window.
The size of the search window is given by the appropriate value in the
\code{width} vector.
The function is intended to count DNA k-mers in selected regions (e.g. exons)
on DNA sequence.
}
\usage{
countDnaKmers(dna,k,start,width)
}
\arguments{
    \item{dna}{\code{character}. Single DNA sequence (vector of length 1). 
        \code{dna} must not contain other characters than "ATCGN". 
        Capitalization does not matter.
        When a 'N' character is found, the current DNA k-mer is skipped.}
    %
    \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs.}
    %
    \item{start}{\code{numeric}. Vector of (1-based) start positions for
        reading frames.
        Reading frame is counted to the right side of the DNA string.}
    %
    \item{width}{\code{numeric}. Defines size of search window for each start
    position. Must have the same length as \code{start} or length 1
    (in which case the values of \code{width} are recycled.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The start positions for counting of DNA k-mers are all positions in
    \{start,...,start+width-1\}. 
    As the identification of a DNA k-mer scans a sequence window of size k,
    the last allowed start position counting a k-mer is \code{nchar(dna)-k+1}.
    The function throws the error \code{'Search region exceeds string end'} 
    when a value  \code{start + width + k > nchar(dna) + 2} occurs.}
%
\value{\code{matrix}. Each colum contains the motif-count values for one frame.
    The column names are the values in the start vector.
    Each row represents one DNA motif.
    The DNA sequence of the DNA motif is given as row.name.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
seq <- "ATAAATA"
countDnaKmers(seq, 2, 1:3, 3)
}
\seealso{countGenomeKmers}
\keyword{countDnaKmers}
