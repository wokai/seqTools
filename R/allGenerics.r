
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  File      :   allGenerics.r                                              ##
##  Content   :   All static variables and (not directly object related )    ##
##                function declarations                                      ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Changing Fastqq object structure
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
setGeneric("mergeFastqq", function(lhs, rhs) standardGeneric("mergeFastqq"))
setGeneric("meltDownK", function(object, newK) standardGeneric("meltDownK"))

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Preparation for Hierarchical clustering (HC)
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
setGeneric("cbDistMatrix", function(object, nReadNorm = max(nReads(object))) 
                  standardGeneric("cbDistMatrix"))

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Slot accessors
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setGeneric("fileNames",   function(object) standardGeneric("fileNames"))
setGeneric("collectTime", function(object) standardGeneric("collectTime"))
setGeneric("collectDur", function(object) standardGeneric("collectDur"))
setGeneric("getK", function(object) standardGeneric("getK"))
setGeneric("nFiles", function(object) standardGeneric("nFiles"))
setGeneric("nNnucs", function(object) standardGeneric("nNnucs"))
setGeneric("nReads", function(object) standardGeneric("nReads"))
setGeneric("maxSeqLen", function(object) standardGeneric("maxSeqLen"))
setGeneric("seqLenCount", function(object) standardGeneric("seqLenCount"))
setGeneric("nucFreq", function(object, i) standardGeneric("nucFreq"))
setGeneric("gcContent", function(object, i) standardGeneric("gcContent"))

setGeneric("gcContentMatrix", function(object)
                                    standardGeneric("gcContentMatrix"))

setGeneric("seqLen", function(object) standardGeneric("seqLen"))
setGeneric("kmerCount", function(object) standardGeneric("kmerCount"))

setGeneric("probeLabel", function(object, label)
                                        standardGeneric("probeLabel"))

setGeneric("probeLabel<-", function(object, value)
                                        standardGeneric("probeLabel<-"))

setGeneric("phred", function(object, i) standardGeneric("phred"))


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Retrieving Phred distribution
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setGeneric("phredQuantiles", function(object, quantiles, i, ...)
                                        standardGeneric("phredQuantiles"))

setGeneric("plotPhredQuant", function(object, i, main, ...)
                                        standardGeneric("plotPhredQuant"))

## Global Phred content functions
setGeneric("phredDist", function(object, i) standardGeneric("phredDist"))

setGeneric("plotPhredDist", function(object, i, maxp=45, col, ...)
                                        standardGeneric("plotPhredDist"))

setGeneric("propPhred", function(object, greater=30, less=93)
                                            standardGeneric("propPhred"))


setGeneric("mergedPhred", function(object) standardGeneric("mergedPhred"))

setGeneric("mergedPhredQuantiles", function(object, quantiles)
                                    standardGeneric("mergedPhredQuantiles"))

setGeneric("plotMergedPhredQuant", function(object, main, ...)
                                    standardGeneric("plotMergedPhredQuant"))

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Predefined Plot functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setGeneric("plotNucFreq", function(object, i, main, maxx, ...)
                                            standardGeneric("plotNucFreq"))

setGeneric("plotGCcontent", function(object,main,...)
                                        standardGeneric("plotGCcontent"))

setGeneric("plotNucCount", function(object, nucs=16, maxx, ...)
                                        standardGeneric("plotNucCount"))

setGeneric("plotKmerCount",
                    function(object, index, mxey, main="K-mer count", ...)
                                            standardGeneric("plotKmerCount"))



## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
