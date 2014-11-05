
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  File      :   fastaFunctions.r                                           ##
##  Content   :   Functions which work on fastq and fastq and write output   ##
##                files: writeFai,  countFastaKmers                          ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

writeFai <- function(infiles, outfiles)
{
    if(!is.character(infiles))
        stop("infiles must be character.")
    
    if(!is.character(outfiles))
        stop("outfiles must be character.")
    
    if(length(infiles) != length(outfiles))
        stop("infiles and outfiles must have equal length.")
    
    if(any( !file.exists(infiles) ) )
        stop("Some infile(s) do not exist!\n")
    
    .Call("write_fai", infiles, outfiles, PACKAGE = "seqTools")
    message("[write_fai] done.\n")
    
    return(invisible())
}


countFastaKmers <- function(filenames, k=4)
{
    if(!is.character(filenames))
        stop("filenames must be character.")
    
    if(!is.numeric(k))
        stop("k must be numeric.")
    k <- as.integer(k)
    
    if( (k < 0) || (k > max_k) )    
        stop("k must be in range 0,    ... , ", max_k, ".")
    
    res <- .Call("count_fasta_Kmers", filenames, k, PACKAGE="seqTools")
    
    return(res)
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
