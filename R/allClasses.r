## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  Version   :   0.9.0                                                      ##
##  File      :   allClasses.r                                               ##
##  Content   :   Declarations for all S4 classes in package                 ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

.onUnload<-function(libpath) { library.dynam.unload("seqTools",libpath) }

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
setClass("Fastqq",representation(
    filenames="character",
    probeLabel="character",
    nFiles="integer",
    k="integer",
    maxSeqLen="integer",
    kmer="matrix",
    firstKmer="matrix",
    nReads="integer",
    nN="integer",
    seqLenCount="matrix",
    gcContent="matrix",
    nac="list",
    phred="list",
    seqLen="matrix",
    collectTime="list"))
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
