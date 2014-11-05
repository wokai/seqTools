\name{Fastqq-class}
\Rdversion{1.1}
\docType{class}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{Fastqq-class}
%
\alias{[-methods}
\alias{[,Fastqq-method}
%
\alias{gcContent}
\alias{gcContent-methods}
\alias{gcContent,Fastqq-method}
%
\alias{getK}
\alias{getK-methods}
\alias{getK,Fastqq-method}
%
\alias{fileNames}
\alias{fileNames-methods}
\alias{fileNames,Fastqq-method}
%
\alias{nFiles}
\alias{nFiles-methods}
\alias{nFiles,Fastqq-method}
%
\alias{nNnucs}
\alias{nNnucs-methods}
\alias{nNnucs,Fastqq-method}
%
\alias{nReads}
\alias{nReads-methods}
\alias{nReads,Fastqq-method}
%
\alias{maxSeqLen}
\alias{maxSeqLen-methods}
\alias{maxSeqLen,Fastqq-method}
%
\alias{phred}
\alias{phred-methods}
\alias{phred,Fastqq-method}
%
\alias{phredQuantiles}
\alias{phredQuantiles-methods}
\alias{phredQuantiles,Fastqq-method}
%
\alias{seqLenCount}
\alias{seqLenCount-methods}
\alias{seqLenCount,Fastqq-method}
%
\alias{nucFreq}
\alias{nucFreq-methods}
\alias{nucFreq,Fastqq-method}
%
\alias{seqLen}
\alias{seqLen-methods}
\alias{seqLen,Fastqq-method}
%
\alias{kmerCount}
\alias{kmerCount-methods}
\alias{kmerCount,Fastqq-method}
% 
\alias{probeLabel}
\alias{probeLabel-methods}
\alias{probeLabel,Fastqq-method}
%
\alias{probeLabel<-}
\alias{probeLabel<--methods}
\alias{probeLabel<-,Fastqq-method}
%
\title{Class \code{"Fastqq"}}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
Contains quality related summarizing data on FASTQ files.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{fastqq("test.fq")}.
}
%
\section{Slots}{
  \describe{
    \item{\code{filenames}:}{\code{"character"}: Vector of Fastqq file names.}
    \item{\code{probeLabel}:}{\code{"character"}: Vector of probe labels.}
    \item{\code{nFiles}:}{\code{"integer"}: Length of fileNamess.}
    \item{\code{k}:}{\code{"integer"}: Length of counted DNA k-mers.}
    \item{\code{maxSeqLen}:}{\code{"integer"} Maximum sequence length found in
        FASTQ files. Determines row-number in 'seqLenCount' matrix and
        column-number in 'nac' and 'phred' slot.}
    \item{\code{kmer}:}{\code{"matrix"} Matrix containing DNA k-mers counts.}
    \item{\code{firstKmer}:}{\code{"matrix"} Matrix containing count of
        incipient DNA k-mers.}
    \item{\code{nReads}:}{\code{"integer"} Vector containing number of reads
        per file.}
    \item{\code{seqLenCount}:}{\code{"matrix"} Matrix containing Counts of
        read lengths.}
    \item{\code{gcContent}:}{\code{"matrix"} Matrix containing GC content
        (in percent).}
    \item{\code{nN}:}{\code{"integer"} Vector containing Number of N 
        nucleotide entries per file.}
    \item{\code{nac}:}{\code{"list"} Contains counted per position alphabet
        frequencies.}
    \item{\code{phred}:}{\code{"list"} Contains per position phred count
        tables (one per Fastqq file).}
    \item{\code{seqLen}:}{\code{"matrix"} Contains minimal and maximal 
        sequence length (one column per file).}
    \item{\code{collectTime}:}{\code{"list"} Contains start and end time of 
        FASTQ reading as 'POSIXct'.}
  }
}
%
\section{Methods}{
The following methods are defined for class \code{Fastqq}:
    \describe{
    Basic accessors:
        \item{getK}{\code{signature(object="Fastqq")}: Returns k-value
            (length of DNA k-mers) as \code{integer}.}
        %
        \item{kmerCount}{\code{signature(object="Fastqq")}: Returns 
            \code{matrix} with 4^k rows anc \code{nFiles} columns. For each
            k-mer and FASTQ-file, the absolute count value of the k-mer in the
            FASTQ file is given.}
        %
        \item{nFiles}{\code{signature(object="Fastqq")}: Returns number of
            Files from which data has been collected as \code{integer}.}
        %
        \item{nNnucs}{\code{signature(object="Fastqq")}: Returns \code{integer}
            vector of length \code{nFiles}. For each FASTQ file, the absolute
            number of containes 'N' nucleotide entries is given.}
        %
        \item{nReads}{\code{signature(object="Fastqq")}: Returns number of 
            reads in each FASTQ file as \code{integer}.}
        %
        \item{fileNames}{\code{signature(object="Fastqq")}: Returns number 
            names of FASTQ files from which data has been collected as 
            \code{character}.}
        %
        \item{maxSeqLen}{\code{signature(object="Fastqq")}: Returns maximum
            sequence length which has been found in all FASTQ files as
            \code{integer}.}
        %
        \item{seqLenCount}{\code{signature(object="Fastqq")}: Returns matrix
            which tables counted read length in all FASTQ files.}
        %
        \item{gcContent}{\code{signature(object="Fastqq",i="numeric")}: 
            Returns \code{integer} vector of length 100 which countains
            absolute read count numbers for each percentage of GC-content.
            \code{i} is the index of the FASTQ file for wich the values
            are returned. The GC content values for all files together can
            be obtained using \code{gcContentMatrix}.}
        %
        \item{nucFreq}{\code{signature(object="Fastqq",i="integer")}:
            Returns \code{matrix} which contains the absolute nucleotide 
            count values for each nucleotide and read position. \code{i} is 
            the index of the FASTQ file for wich the values are returned.}
        %
        \item{seqLen}{\code{signature(object="Fastqq")}: Returns \code{matrix}
            with two rows and \code{nFiles} columns. For each file the minimum
            and maximum read length is given.}
        %
        \item{kmerCount}{\code{signature(object="Fastqq")}: Returns a 
            \code{matrix} with 4^k rows and \code{nFiles} columns. Each entry 
            gives the absolute count of the k-mer (given as row name) in each
            file (given as column name).}
        %
        \item{phred}{\code{signature(object="Fastqq",i="integer")}: Returns
            a \code{matrix} with 93 rows and \code{maxSeqLen} columns. The 
            matrix gived the absolute counts of each phred value for each 
            sequence position. \code{i} is the index of the FASTQ file for 
            wich the values are returned.}
        %
        \item{phredQuantiles}{\code{signature(object="Fastqq", 
            quantiles="numeric", i="integer")}: Returns a \code{data.frame}.
            The data.frame has one row for each given quantile and
            \code{maxSeqLen} columns. Each value gives the quantile (given by
            row name) of the phred values at the sequence position (given by
            column name). For the \code{quantiles} argument, a numeric vector 
            with values in [0,1] must be given. For the \code{i} argument, a 
            single integer value must be given which denotes the index of the
            FASTQ file from which values are returned (value must be in
            \{1,...,nFiles\}).}
        %
        \item{probeLabel}{\code{signature(object="Fastqq")}: Returns 
            \code{character} vector which contains the \code{probeLabel} 
            entries for given \code{Fastqq} object.}
    }
}
%
\references{
Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM
The sanger FASTQ file format for sequences with quality scores and
    the Solexa/Illumina FASTQ variants.
Nucleic Acids Research 2010 Vol.38 No.6 1767-1771
}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
basedir <- system.file("extdata", package="seqTools")
setwd(basedir)
fq <- fastqq(c("g4_l101_n100.fq.gz","g5_l101_n100.fq.gz"), 
                                k=4, probeLabel=c("g4","g5"))
#
fileNames(fq)
getK(fq)
nNnucs(fq)
nFiles(fq)
nReads(fq)
maxSeqLen(fq)
collectTime(fq)
collectDur(fq)
slc<-seqLenCount(fq)
nf<-nucFreq(fq,1)
nf[1:4,1:10]
seqLen(fq)
probeLabel(fq)
probeLabel(fq) <- 1:nFiles(fq)
#
kc<-kmerCount(fq)
kc[1:10, ]
plotKmerCount(fq)
#
ph<-phred(fq, 1)
ph[25:35,1:15]
pq <- phredQuantiles(fq,c(0.25, 0.5, 0.75), 1)
plotNucFreq(fq, 1)
# Nucleotide count
plotNucCount(fq, 2:3) 
# GC content
gcContent(fq, 1)
#
fqq<-fq[1]
}
\seealso{fastqq}
\keyword{classes}
\keyword{fastqq}
\keyword{kmer}
