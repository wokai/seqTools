\name{countSpliceKmers}
%
\alias{countSpliceKmers}
%
\title{countSpliceKmers: Counting K-mers on donor (5', upstream) sides
    (exonic) of splice sites.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
The function regards the given string as DNA sequence bearing a collection of
    splice sites. The given \code{lEnd} and \code{rStart} positions act as 
    (1-based) coordinates of the innermost exonic nucleotides. They reside on
    exon-intron boundaries and have one exonic and one intronic adjacent
    nucleotide. The function counts \code{width} k-mers upstream on exonic
    DNA in reading direction (left -> right on (+) strand, right -> left on
    (-) strand).
}
%
\usage{
countSpliceKmers(dna, seqid, lEnd, rStart, width, strand, k)
}
\arguments{
    \item{dna}{\code{character}. Vector of DNA sequences. \code{dna} must
        not contain other characters than "ATCGN".
        Capitalization does not matter.
        When a 'N' character is found, the current DNA k-mer is skipped.}
    %
    \item{seqid}{\code{numeric}. Vector of (1-based) values coding for one
        of the given sequences.}
    %
    \item{lEnd}{\code{numeric}. Vector of (1-based) left-end positions.
        Will be used as rightmost window position.}
    %
    \item{rStart}{\code{numeric}. Vector of (1-based) right-start positions.
        Will be used as leftmost window positions (over which(n-1) positions
        overhang will be counted as part of frame).}
    %
    \item{width}{\code{numeric}. Vector of window width values.}
    %
    \item{strand}{\code{factor or numeric}. First factor level (or numeric: 1)
        value will be interpreted as (+) strand
        For any other values, the reversed complement sequence will be
        counted (in left direction from start value).
        For (+) strand, the lEnd value will be used as starting position.
        For (-) strand, the rStart position will be used as starting
        positions.}
    %
    \item{k}{\code{numeric}. Number of nucleotides in tabled DNA motifs.
        Only a single value is allowed (length(n) = 1 !)}  
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
seq <- "acgtGTccccAGcccc"
countSpliceKmers(seq, seqid=1, lEnd=4, rStart=10, width=2, strand=1, k=3)
#
sq1 <- "TTTTTCCCCGGGGAAAA"
sq2 <- "TTTTTTTCCCCGGGGAAAA"
sq <- c(sq1, sq2)
seqid <- c( 1, 1, 2, 2)
lEnd  <- c( 9, 9, 11, 11)
rStart <- c(14, 14, 16, 16)
width <- c( 4, 4, 4, 4)
strand <- c( 1, 0, 1, 0)
countSpliceKmers(sq, seqid, lEnd, rStart, width, strand, k=2)
}
\keyword{countSpliceKmers}
