
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  File      :   trimFastq.r                                                ##
##  Content   :   Functions which work on fastq and fastq and write output   ##
##                files: trimFastq, writeFai                                 ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## trimFastq: Trimming and discarding reads based on quality values
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

trimFastq<-function(infile, outfile="keep.fq.gz", discard="disc.fq.gz", 
        qualDiscard=0, qualMask=0, fixTrimLeft=0, fixTrimRight=0,
        qualTrimLeft=0, qualTrimRight=0, qualMaskValue=78, minSeqLen=0)
{
    if(!is.character(infile))
        stop("'infile' must be a string.")
    
    if(length(infile) != 1)
        stop("'infile' must have length 1.")
    
    if(!file.exists(infile))
        stop("File 'infile' not found.")
    
    if(!is.character(outfile))
        stop("'outfile' must be character.")
    
    if(!is.character(outfile))
        stop("'discard' (file) must be character.")
    
    if(!is.character(discard))
        stop("'discard' must be character.")
    
    if(!is.character(discard))
        stop("'discard' (file) must be character.")
    
    if(!is.numeric(qualDiscard))
        stop("'qualDiscard' must be numeric.")
    
    if(length(qualDiscard)!=1)
        stop("'qualDiscard' must have length 1.")
    
    qualDiscard<-as.integer(qualDiscard)
    if( (qualDiscard < 0) || (qualDiscard > 93))
        stop("'qualDiscard' out of range.")
    
    if(!is.numeric(qualMask))
        stop("'qualMask' must be numeric.")
    
    if(length(qualMask) != 1)
        stop("'qualMask' must have length 1.")
    
    qualMask<-as.integer(qualMask)
    if( (qualMask < 0) || (qualMask > 93) )
        stop("'qualMask' out of range.")
    
    if(!is.numeric(fixTrimLeft))
        stop("'fixTrimLeft' must be numeric.")
    
    if(length(fixTrimLeft) != 1)
        stop("'fixTrimLeft' must have length 1.")
    
    if( (fixTrimLeft < 0) || (fixTrimLeft > 100) )
        stop("'fixTrimLeft' out of range.")
    
    if(!is.numeric(fixTrimRight))
        stop("'fixTrimRight' must be numeric.")
    
    if(length(fixTrimRight) != 1)
        stop("'fixTrimRight' must have length 1.")
    
    if( (fixTrimRight < 0) || (fixTrimRight > 100) )
        stop("'fixTrimRight' out of range.")
    
    if(!is.numeric(qualTrimLeft))
        stop("'qualTrimLeft' must be numeric.")
    
    if(length(qualTrimLeft)!= 1)
        stop("'qualTrimLeft' must have length 1.")
    qualTrimLeft<-as.integer(qualTrimLeft)
    
    if( (qualTrimLeft < 0) || (qualTrimLeft > 93) )
        stop("'qualTrimLeft' out of range.")
    
    if(!is.numeric(qualTrimRight))
        stop("'qualTrimRight' must be numeric.")
    
    if(length(qualTrimRight) != 1)
        stop("'qualTrimRight' must have length 1.")
    qualTrimRight<-as.integer(qualTrimRight)
    
    if( (qualTrimRight < 0) || (qualTrimRight > 93) )
        stop("'qualTrimRight' out of range.")
    
    if(!is.numeric(qualMaskValue))
        stop("'qualMaskValue' must be numeric.")
    
    if(length(qualMaskValue) != 1)
        stop("'qualMaskValue' must have length 1.")
    
    if( (qualMaskValue < 0) || (qualMaskValue > 93) )
        stop("'qualMaskValue' out of range.")
    
    if(!is.numeric(minSeqLen))
        stop("'minSeqLen' must be numeric.")
    
    if(length(minSeqLen) != 1)
        stop("'minSeqLen' must have length 1.")
    
    if( (minSeqLen < 0) || (minSeqLen > 200) )
        stop("'minSeqLen' out of range.")
    
    val<-as.integer(c(
        fixTrimLeft,
        fixTrimRight,
        qualTrimLeft,
        qualTrimRight,
        qualDiscard,
        qualMask,
        qualMaskValue,
        minSeqLen
    ))
    res<-.Call("trim_fastq", infile, val, c(outfile,discard))
    
    bm <- Sys.localeconv()[7]
    message("[trimFastq] ", format(res[1], width=11, big.mark=bm),
                                            " records written to outfile.")
    message("[trimFastq] ", format(res[2], width=11, big.mark=bm),
                                            " records written to discard.")
    return(invisible(res))
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
