\name{revCountDnaKmers}
%
\alias{revCountDnaKmers}
%
\title{revCountDnaKmers: Counting K-mers in DNA sequence.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Desription
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Counts DNA K-mers for reverse complement of given DNA sequence.
    The k-mers are counted in a set of search windows,
    which are defined by \code{start} and \code{width} parameter.
    From each position of the search window, a DNA k-mer is identified on the
    left hand side on the reverse complement of the given DNA sequence.
    Each value in the \code{start} vector defines the right border
    of a search window.
    The size of the search window is given by the appropriate value
    in the \code{width} vector. 
}
%
\usage{revCountDnaKmers(dna, k ,start, width)}
%
\arguments{
    \item{dna}{\code{character}. Single DNA sequence (vector of length 1).
        \code{dna} must not contain other characters than "ATCGN".
        Capitalization does not matter.
        When a 'N' character is found, the ongoing identification
        of a DNA k-mer is terminated.}
    %
    \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs.}
    %
    \item{start}{\code{numeric}. Vector of (1-based) start positions for 
        reading frames.}
    %
    \item{width}{\code{numeric}. Defines number of k-mers (size of search
        window) for each start position. Must have the same length as 
        \code{start} or length 1 (in which case the values of \code{width} 
        are recycled.)}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The start positions for identification of DNA k-mers are all positions
    in \{start-width+1,...,start\}. In order to prevent counting before the
    first nucleotide of the DNA sequence, all start values must be >= width + k.
    The function throws an error when this border is exceeded.
}
%
\value{\code{matrix}. Each colum contains the motif-count values for one frame.
    Each row represents one DNA motif. The DNA sequence of the DNA motif is
    given as row.name.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
rseq <- "TATTAT"
revCountDnaKmers(rseq, 2,6:4, 2)
}
%
\keyword{revCountDnaKmers}
