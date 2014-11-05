
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  File      :   allStatics.r                                               ##
##  Content   :   All static variables and (not directly object related )    ##
##                function declarations                                      ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Global variables:
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

strandlevels <- c("+", "-", "*")
iupac_chars <- c('a', 'c', 'g', 't', 'u', 'r', 'y', 's', 'w', 'k', 'm', 'b',
                'd', 'h', 'v', 'n', '.', '-', '=')

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Although static array sizes would allow larger k (e.g. 15: 8 GB RAM), 
## exponential increase in run-time restricts usability to k values in 1:12
## 4^12 = 16.777.216 possible combinations.
## 'Hard' coded (in stat_defs.h: 15)
max_k <- 12
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## End global variables
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Some useful functions for work with phred's
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

char2ascii <- function(c) { strtoi(charToRaw(c), base=16L) }

ascii2char <- function(x, multiple=FALSE)
                                { rawToChar(as.raw(x), multiple=multiple) }

phredTable <- function(phred=0:93)
{
    if(!is.numeric(phred))
        stop("phred must be numeric.")
    
    phred <- sort(unique(as.integer(phred)))
    
    if( (phred[1] < 0) || (max(phred) > 93) )
        stop("Only phred values in 0:93 are allowed.")

    ascii <- phred + 33
    return(data.frame(ascii=ascii, phred=phred,
                    char=ascii2char(ascii, multiple=TRUE)))
}
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

# unexported functions,  do not check for x
rel_int  <- function(x) {return(.Call("rel_int", x, PACKAGE="seqTools"))}
rel_real <- function(x) {return(.Call("rel_real", x, PACKAGE="seqTools"))}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
enlargeIntMat <- function(x, newDim){
    if(!is.matrix(x))
        stop("x must be matrix.")

    if(!is.numeric(newDim))
        stop("newDim must be numeric.")
    
    newDim <- as.integer(newDim)

    if(length(newDim) != 2)
        stop("newDim must have length 2.")

    return(.Call("r_enlarge_int_mat", x, newDim, PACKAGE="seqTools"))
}

kMerIndex <- function(kMers, k=nchar(kMers)[1], base=1)
{
    if(!is.character(kMers))
        stop("kMers must be character.")

    if(!is.numeric(k))
        stop("k must be numeric.")
    k <- as.integer(k)

    if(k<0 || k > max_k)
        stop("k must be in range 0, ... , ", max_k, ".")

    if(!all(nchar(kMers) == k))
        stop("All kMers must have k characters!")

    if(!is.numeric(base))
        stop("base must be numeric")

    if(length(base) > 1)
        stop("base must have length 1.")
    base<-as.integer(base)

    if(!(base %in% 0:1))
        stop("base must be 0 or 1.")

    return(.Call("get_Kmer_Index", kMers, k, PACKAGE="seqTools") + base)
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
