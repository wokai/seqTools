
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  File      :   kMer.r                                                     ##
##  Content   :   Functionality for counting DNA k-mers                      ##
##                (independent of fastq or fasta files)                      ##
##                countKmers, countDnaKmers, revCountDnaKmers,               ##
##                countGenomeKmers, countSpliceKmers                         ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on specified regions inside single (character) sequence
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
countDnaKmers <- function(dna, k, start=1, width=nchar(dna) - k + 1)
{
    if(!is.character(dna))
        stop("'dna' must be character.")
    
    if(length(dna) != 1)
        stop("'dna' must have length 1.")
    
    if(is.numeric(start))
        start <- as.integer(start)
    
    if(is.numeric(width))
        width <- as.integer(width)
    
    if(length(width) == 1)
        width <- rep(width, length(start))
    
    if(is.numeric(k))
        k <- as.integer(k)
    
    if(length(k)!=1)
        stop("'k' must have length 1.")
    
    if(k < 1)
        stop("'k' must be positive.")
    
    if(k > max_k)
        stop("'k' must not exceed", max_k, ".")
  
    nc <- nchar(dna)
    if(k > nc)
        stop("'k' must be <= nchar(dna).")
    if(any(start + width + k > nc + 2))
        stop("Search region exceeds string end.")
  
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Counts N's
    ## ToDo: Return value
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    nn <- integer(length(start))
    return(.Call("count_dna_Kmers", 
                            dna, start, width, k, nn, PACKAGE="seqTools"))
    
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on specified regions inside single (character) sequence 
## in reverse direction
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
revCountDnaKmers <- function(dna, k, start, width)
{
    if(!is.character(dna))
        stop("'dna' must be character.")
    
    if(length(dna) != 1)
        stop("'dna' must have length 1.")
    
    if(is.numeric(start))
        start <- as.integer(start)
    
    if(is.numeric(width))
        width <- as.integer(width)
    
    if(length(width) == 1)
        width <- rep(width, length(start))
    
    if(length(width) != length(start))
        stop("'width' must have length 1 or the same length as 'start'.")
    
    if(is.numeric(k))
        k <- as.integer(k)
    
    if(any(width + k > start))
        stop("'width' must be <=  'start' - 'k'.")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## nn contains N counts
    ## ToDo: Add value of nn to returned object
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    nn <- integer(length(start))
    
    return(.Call("rev_count_dna_Kmers",
                            dna, start, width, k, nn, PACKAGE="seqTools"))
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on specified regions inside multiple (character) sequences
## in possibly reversed direction (depending on strand)
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
countGenomeKmers <- function(dna, seqid, start, width, strand, k)
{
    if(!is.character(dna))
        stop("'dna' must be character.")
    
    if(!is.numeric(seqid))
        seqid <- as.integer(seqid)
    rg <- range(seqid)
    
    if(rg[1] < 0)
        stop("Negative seqid's are not allowed.")
    
    if(rg[2] > length(dna))
        stop("Out of range seqid's.")
    
    if(!is.numeric(start))
        stop("'start' must be numeric.")
    start <- as.integer(start)
    
    if(!is.numeric(width))
        stop("'width' must be numeric")
    width <- as.integer(width)
    
    if(is.factor(strand))
        strand <- as.integer(strand)
    else
    {
        if(!is.numeric("strand"))
            strand <- as.integer(strand)        
    }
    
    nStart <- length(start)
    if( (length(seqid) != nStart) | (length(width) != nStart) |
                (length(strand) != nStart) )
        stop("'seqid', 'start', 'width' and 'strand' must have same length.")
    
    if(length(k) != 1)
        stop("'k' must be a single value.")
    
    if(!is.numeric(k))
        stop("'k' must be numeric.")
    k<-as.integer(k)
    
    if(k <= 0)
        stop("'k' must be >=1")
    
    if(k > max_k)
        stop("'k' must not exceed", max_k, ".")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Counts N's
    ## ToDo: Return value
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    nn <- integer(length(start))
    return(.Call("count_genome_Kmers", dna, seqid, start, width,
                            strand, k, nn, PACKAGE="seqTools"))
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Counts DNA k-mers on each border of a splice-site defined by wLend and 
## wRstart in range of size width
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

countSpliceKmers <- function(dna, seqid, lEnd, rStart, width, strand, k)
{
    if(!is.character(dna))
        stop("'dna' must be character.")
    
    if(!is.numeric(seqid))
        stop("'seqid' must be numeric.")
    seqid <- as.integer(seqid)
    
    rg <- range(seqid)
    if(rg[1] < 0)
        stop("Negative seqid's are not allowed.")
    
    if(rg[2] > length(dna))
        stop("Out of range seqid's.")
    
    if(!is.numeric(lEnd))
        stop("'lEnd' must be numeric.")
    lEnd <- as.integer(lEnd)
    
    if(!is.numeric(rStart))
        stop("'rStart' must be numeric.")
    rStart <- as.integer(rStart)
    
    if(!is.numeric(width))
        stop("'width' must be numeric.")
    width <- as.integer(width)
    
    if(is.factor(strand))
        strand <- as.integer(strand)
    else
    {
        if(!is.numeric("strand"))
            strand <- as.integer(strand)        
    }
    
    nStart <- length(lEnd)
    if(length(seqid) != nStart | length(rStart) != nStart | 
             length(width) != nStart | length(strand) != nStart)
    {
        stop(paste("'seqid', 'lEnd', 'rStart', 'width'",
                            "and 'strand' must have equal length."))
    }
    if(!is.numeric(k))
        stop("'k' must be numeric.")
    k <- as.integer(k)
    
    if(length(k) != 1)
        stop("'k' must be a single value.")
    
    if(k > max_k)
        stop("'k' must not exceed", max_k, ".")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Plus strand
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    plus_strand<-strand == 1
    
    if(sum(plus_strand) > 0)
    {
        if(any((lEnd[plus_strand] - width[plus_strand] - k + 1) < 0))
            stop("lEnd must be >= width+k-1 for all +-strand coordinates")
    }
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Counts N's
    ## ToDo: Return value
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    nn <- integer(length = nStart)
    
    return(.Call("count_splice_Kmers", dna, seqid, lEnd, rStart, width, 
                            strand, k, nn, PACKAGE="seqTools")) 
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
