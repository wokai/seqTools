\name{cbDistMatrix}
%
\alias{cbDistMatrix}
\alias{cbDistMatrix-methods}
\alias{cbDistMatrix,Fastqq-method}
%
\title{
cbDistMatrix function: Calculates pairwise distance matrix from DNA k-mer counts
based on a modified Canberra distance.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
Calculates pairwise distance matrix from DNA k-mer counts based on a modified
Canberra distance. Before calculating canberra distances, read counts are 
normalized (in order to correct systematic effects on the distance) by 
scaling up read counts in each DNA k-mer count vector so that normalized
read counts in each sample are nearly equal.
}
\usage{
cbDistMatrix(object,nReadNorm=max(nReads(object)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{\code{Fastqq}: Object from which DNA k-mer counts are used.}
    \item{nReadNorm}{\code{numeric}: 
    Number of reads per file to wich all contained DNA k-mer counts are
    normalized. 
    Because the normalization is intended to increase counts the value
    must be greater than all FASTQ file read counts (as reported by nReads).
    Therefore the standard value is chosen to the maximal number of reads
    recorded in this object.
    This normalization is necessary to compensate for systematic effects 
    in the canberra distance.}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The distance between two DNA k-mer normalized count vectors
        is calculated by
            \deqn{df (X,Y)= \sum_{i=1}^n cbd(x_i,y_i) / 4^k}{%
                df (X,Y) = \sum cbc(x_i, y_i) / 4^k}
        where cb is given by
            \deqn{cbd(x,y)=|x-y|/(x+y).}
}
%
\value{Square \code{matrix}. The number of rows equals the number of files 
                (=\code{nFiles(object)}). }
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and the
Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
\author{
Wolfgang Kaisers
}
\note{The static size of the retured k-mer array is 4^k.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir<-system.file("extdata",package="seqTools")
basenames<-c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz")
filenames<-file.path(basedir,basenames)
fq<-fastqq(filenames,6,c("g4","g5"))
dm<-cbDistMatrix(fq)
}
\keyword{cbDistMatrix}
\keyword{kmer}
\seealso{hclust}
