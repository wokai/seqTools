
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  Content   :   Doing some diagnostic and interventional tasks on fastq    ##
##                and fasta                                                  ##
##                esp. DNA k-mer counts.                                     ##
##  Version   :   0.99.34                                                    ##
##                                                                           ##
##  Changelog :                                                              ##
##  26.08.13  :   0.0.1 Project created                                      ##
##  03.09.13  :   0.0.6 C-Code valgrind tested                               ##
##  27.09.13  :   0.1.0 Added fastqDnaKmers                                  ##
##  14.10.13  :   0.1.1 Added C-library for parsing gz fasta and fastq files ##
##  17.10.13  :   0.1.3 C-Wrapper for fastq working.                         ##
##  17.10.13  :   0.1.6 First version of R package                           ##
##  21.10.13  :   0.3.0 New C library for fastq parsing                      ##
##  28.10.13  :   0.4.0 Added fastq-loc functions                            ##
##  29.10.13  :   0.4.4 seqTools C-code valgrind tested.                     ##
##  01.11.13  :   0.5.0 Distance matrices implemented                        ##
##  02.11.13  :   0.5.1 First working version with clustering based on       ##
##                      K-mers                                               ##
##  07.11.13  :   0.5.4 countFastaKmers now resizes arrays automatically     ##
##  09.11.13  :   0.6.0 count_fastq now resizes arrays automatically         ##
##  11.11.13  :   0.6.5 Added fastq simulation functions                     ##
##  19.11.13  :   0.7.0 Added trimFastq function                             ##
##  30.11.13  :   0.9.0 Separated R source files                             ##
##  22.12.13  :   0.99.2 Added '['-operator for Fastqq class                 ##
##  19.01.14  :   0.99.7 Added zlib version control for correction of        ##
##                       gzBuffer                                            ##
##                       error                                               ##
##                        + checks: cran win-builder + valgrind              ##
##  21.05.14  :   0.99.34 Corrected error in count_fasta_Kmers               ##
##                        which freezed function                             ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

.onUnload <- function(libpath) { library.dynam.unload("seqTools", libpath) }

## see: http://bioconductor.org/developers/how-to/coding-style/

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Data collection functions:
## Fastqq,  fastqKmerLocs,  fastqKmerSubsetLocs
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

fastqq <- function(filenames, k=6, probeLabel)
{
    k <- as.integer(k)
    
    tl <- list()
    tl$start <- Sys.time()
    filenames <- path.expand(filenames)
    
    res <- .Call("count_fastq", filenames, k, PACKAGE="seqTools")
    
    tl$end <- Sys.time()
    res@collectTime <- tl
    
    # Correct minSeqLen when empty files are counted
    if(any(res@nReads==0))
        res@seqLen[1,res@nReads==0] <- 0
    
    if(!missing(probeLabel))
    {
        if(length(probeLabel) == res@nFiles)
            res@probeLabel <- as.character(probeLabel) 
        else{
            warning("[Fastqq] probeLabel and filenames must have equal length.")
            res@probeLabel <- as.character(1:res@nFiles)
        }
    }else{
        res@probeLabel <- as.character(1:res@nFiles)
    }
    return(res)
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## fastq K-mer locs
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

fastqKmerLocs <- function(filenames, k=4)
{
    if(!is.numeric(k))
        stop("'k' must be numeric.")
    k <- as.integer(k)
    
    if( (k < 0) || (k > max_k) )
        stop("'k' must be in range 0, ... , 16.")
    
    return(.Call("fastq_Kmer_locs", filenames, k, PACKAGE="seqTools"))
}


fastqKmerSubsetLocs <- function(filenames, k=4, kIndex)
{
    # Returns list with matrix elements.
    if(!is.numeric(k))
        stop("'k' must be numeric.")
    
    k <- as.integer(k)
    
    if( (k < 0) || (k > max_k) )
        stop("'k' must be in range 0, ... , ", max_k, ".")
    
    if(!is.numeric(kIndex))
        stop("'kIndex' must be numeric.")
    kIndex <- as.integer(kIndex)
    
    if(any(kIndex < 0))
        stop("No negative 'kIndex' values allowed.")
    
    if(any(kIndex > (4^k)) )
        stop("'kIndex' out of range (>4^k).")
    
    return(.Call("fastq_KmerSubset_locs",
                                filenames, k, kIndex, PACKAGE="seqTools"))
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## End: Data collection functions.
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Standard slot accessor functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


setMethod("fileNames", "Fastqq", function(object)
                                            {return(object@filenames)})

setMethod("collectTime", "Fastqq", function(object)
                                            {return(object@collectTime)})

setMethod("collectDur", "Fastqq", function(object) {
return(as.numeric(difftime(object@collectTime$end, object@collectTime$start,
                                                    units = "secs")))
})

setMethod("getK", "Fastqq", function(object) {return(object@k)})

setMethod("nFiles", "Fastqq", function(object) {return(object@nFiles)})

setMethod("nNnucs", "Fastqq", function(object) {return(object@nN)})

setMethod("nReads", "Fastqq", function(object) {return(object@nReads)})

setMethod("maxSeqLen", "Fastqq", function(object) {return(object@maxSeqLen)})

setMethod("seqLenCount", "Fastqq", function(object)
{
    res<-object@seqLenCount
    colnames(res) <- object@probeLabel
    return(res)
})

setMethod("nucFreq", "Fastqq", function(object, i)
{
    if(missing(i))
        stop("Argument 'i' is not optional.")
    
    if(!is.numeric(i))
        stop("'i' must be numeric.")
    
    i <- as.integer(i)
    if( (i < 1) || (i > object@nFiles) )
        stop("'i' must be >0 and < nFiles.")
    
    return(object@nac[[i]])
})

setMethod("gcContent", "Fastqq", function(object, i)
{
    if(missing(i))
        stop("Argument 'i' is not optional.")
    
    if(!is.numeric(i))
        stop("'i' must be numeric.")
    
    i <- as.integer(i)
    if( (i < 1) || (i > object@nFiles))
        stop("'i' must be >0 and < nFiles.")
    
    return(object@gcContent[, i])
})


setMethod("gcContentMatrix", "Fastqq", function(object)
{
    gcc <- object@gcContent
    colnames(gcc) <- object@probeLabel
    return(gcc)
})

setMethod("seqLen", "Fastqq", function(object)
{
    sql <- object@seqLen
    colnames(sql) <- object@probeLabel
    return(sql)
})

setMethod("kmerCount", "Fastqq", function(object)
{
    kmer <- object@kmer
    colnames(kmer) <- object@probeLabel
    return(kmer)
})


setMethod("probeLabel", "Fastqq", function(object){return(object@probeLabel)})
setReplaceMethod("probeLabel", "Fastqq", function(object, value)
{
    if(length(value) != nFiles(object))
        stop("'value' must have length ", nFiles(object), ".")
    
    val <- as.character(value)
    if(any(table(val)) > 1)
    {
        warning("[probeLabel <- .Fastqq] probeLabel unique suffix added.")
        val <- paste(1:nFiles(object), val, sep="_")
    }
    object@probeLabel <- val
    
    return(object)
})

setMethod("phred", signature="Fastqq", definition=function(object, i)
{
    if(missing(i))
        stop("Argument 'i' is not optional.")
    
    if(!is.numeric(i))
        stop("'i' must be numeric.")
    
    i <- as.integer(i)
    if( (i < 1) || (i > object@nFiles) )
        stop("'i' must be >0 and < nFiles.")
    
    return(object@phred[[i]])
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## End: Standard slot accessor functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


setMethod("show", "Fastqq", function(object)
{
    bm <- Sys.localeconv()[7]
    w <- 20
    r <- "right"
    cat("Class       : ", format(class(object), w=w, j=r)                              , "\n", sep="")
    cat("nFiles      : ", format(format(nFiles(object)          , big.m=bm), w=w, j=r) , "\n", sep="")
    cat("maxSeqLen   : ", format(format(maxSeqLen(object)       , big.m=bm), w=w, j=r) , "\n", sep="")
    cat("k (Kmer len): ", format(format(getK(object)            , big.m=bm), w=w, j=r) , "\n", sep="")
    cat("\n")
    cat("nReads      : ", format(format(sum(as.numeric(nReads(object)))    , big.m=bm), w=w, j=r)   , "\n", sep="")
    cat("nr  N   nuc : ", format(format(sum(nNnucs(object))     , big.m=bm), w=w, j=r) , "\n", sep="")
    cat("Min seq len : ", format(format(min(seqLen(object)[1, ]), big.m=bm), w=w, j=r) , "\n", sep="")
    cat("Max seq len : ", format(format(max(seqLen(object)[2, ]), big.m=bm), w=w, j=r) , "\n", sep="")
    return(invisible())
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Phred related functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("phredQuantiles", "Fastqq", function(object, quantiles, i, ...)
{
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Checking arguments
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    
    ## Check quantiles argument
    if(missing(quantiles))
        stop("'quantiles' argument is not optional")
    
    if(!is.numeric(quantiles))
        stop("Quantiles must be numeric.")
    
    if(!(all(quantiles >= 0) & all(quantiles <= 1) ))
        stop("All quantiles mustbe in [0, 1]")
    
    quantiles <- sort(unique(round(quantiles, 2)))
    
    ## Check 'i' argument
    if(missing(i))
        stop("'i' argument is not optional")
    
    if(!is.numeric(i))
        stop("'i' must be numeric.")
    
    if(length(i) > 1)
        stop("'i' must have length 1")
    
    i <- as.integer(i)
    if( (i < 1) || (i > object@nFiles) )
            stop("'i' must be >0 and <nFiles.")

    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Count qual values for each sequence position
    # Convert integer counts into column-wise relative values
    # Maximum counted read length
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    vec <- 1:seqLen(object)[2, i]
    qrel <- as.data.frame(apply(object@phred[[i]][, vec], 2, rel_int))
    names(qrel) <- vec
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Walk through each column and extract row number
    # for given quantile values
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    res <- .Call("get_column_quantiles", quantiles, qrel, PACKAGE="seqTools")
    return(res)
})


setMethod("plotPhredQuant", "Fastqq", function(object, i, main, ...){
    if(!is.numeric(i))
        stop("'i' must be numeric.")
    
    i <- as.integer(i)
    
    if( (i < 1) || (i > object@nFiles) )
        stop("'i' must be >0 and <nFiles.")
    
    quant <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    cols <- c("#1F78B4", "#FF7F00", "#E31A1C", "#FF7F00", "#1F78B4") 
    qq <- phredQuantiles(object, quant, i)
    maxQ = floor(1.2*max(qq))
    xv <- 1:ncol(qq)
    
    if(missing(main))
    {
        main <- paste("Position wise Phred Quantiles (", 
                                probeLabel(object)[i], ")", sep="")
    }

    plot(xv, xv, ylim=c(0, maxQ), type="n", bty="n", las=1,
        ylab = "Phred score", xlab="Sequence position", main=main, ...)
     
    lines(xv, qq[1, ], col=cols[1], lty=2)
    lines(xv, qq[2, ], col=cols[2], lty=1)
    lines(xv, qq[3, ], col=cols[3], lwd=2)
    lines(xv, qq[4, ], col=cols[4], lty=1)
    lines(xv, qq[5, ], col=cols[5], lty=2)
    
    legend("top", ncol=6, lty=c(2, 1, 1, 1, 2),
            lwd=c(1, 1, 2, 1, 1), col=cols, xjust=0.5,
            legend=c("10%", "25%", "50%", "75%", "90%"), bty="n", cex=0.8)
    return(invisible())
})




## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Global Phred content functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
setMethod("phredDist", "Fastqq", function(object, i){
    idx <- 1:nFiles(object)
    
    if(missing(i))
        i <- idx
    else
    {
        if(!is.numeric(i))
            stop("'i' must be numeric.")
        
        if(!all(is.element(i, idx)))
            stop("'i' is out of range.")
    }
    
    phred <- Reduce("+", object@phred[i])
    phred <- matrix(as.numeric(phred), nrow=nrow(phred))
    
    phred_vals <- apply(phred, 1, sum)
    phred_dist <- phred_vals/sum(phred_vals)
    names(phred_dist) <- rownames(object@phred[[1]])
    
    return(phred_dist)    
})


setMethod("plotPhredDist", "Fastqq", function(object, i, maxp=45, col, ...)
{
    if(!is.numeric(maxp))
        stop("maxp must be numeric")
    
    if(!is.integer(maxp))
        maxp<-as.integer(maxp)
    
    if(maxp <= 0)
        stop("maxp must be >=0")
    
    if(missing(col))
        col <- topo.colors(10)[3]
    
    phred <- phredDist(object, i)
    maxy <- ceiling(max(phred) * 5) / 5
    x <- 1:maxp
    xmax <- 10 * (ceiling(maxp / 10))
    
    plot(x, phred[x], ylim=c(0, maxy), xlim=c(0, xmax), type="l", lwd=2,
            col=col, yaxt="n", bty="n", xlab="Phred value",
            ylab="Content (%)", ...)
    
    ylb <- 0:(10 * maxy) / 10
    axis(side=2, at=ylb, labels=100 * ylb, las=1)
    return(invisible())
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# Not exported:
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
pd_l10 <- function(x){ x <- phredDist(x); return(sum(x[1:10]) / sum(x))}


setMethod("propPhred","Fastqq",function(object, greater=30, less=93)
{
    if(!is.numeric(greater))
        stop("'greater' must be numeric.")
    
    if(length(greater) != 1)
        stop("'greater' must have length 1.")
    
    if(!is.numeric(less))
        stop("'less' must be numeric.")
    
    if(length(less) != 1)
        stop("'less must have length 1.")
    
    ## + + + + + + + + + + + + + + + + + + + + + + ##
    ## greater and less shall be used as
    ## array indices: increase greater
    ## + + + + + + + + + + + + + + + + + + + + + + ##
    greater<-as.integer(greater+1)
    less<-as.integer(less)
    
    if(greater < 1)
        stop("'greater' must be >=0.")
    
    if(less > 93)
        stop("'less' must be < 94.")
    
    if(greater >= less)
        stop("'greater' must be <= 'less'")
    
    n <- nFiles(object)
    res <- numeric(n)
    for(i in 1:n)
    {
        pd <- phredDist(object, i)
        res[i] <- sum(pd[greater:less])
    }
    names(res) <- probeLabel(object)
    return(res)
})

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## End: Global Phred content functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("mergedPhred", "Fastqq", function(object){
    sql <- seqLen(object)
    maxSeqLen <- max(sql[2, ])
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## as.numeric: Sum on integer is likely to exceed 
    ##                 maximal 32-bit integer  values
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    if(sql[2, 1] < maxSeqLen)
    {
        mp <- as.numeric(.Call("r_enlarge_int_mat", object@phred[[1]], 
                c(nrow(object@phred[[1]]), maxSeqLen), PACKAGE="seqTools"))
    }else{
        mp <- as.numeric(object@phred[[1]])
    }
    
    
    n <- nFiles(object)
    if(n > 1)
    {
        for(i in 2:n)
        {
            if(sql[2, i] < maxSeqLen)  
            {
                mp <- mp + as.numeric(.Call("r_enlarge_int_mat",
                            object@phred[[i]],
                            c(nrow(object@phred[[i]]), maxSeqLen), 
                            PACKAGE="seqTools"))
            }else{
                mp <- mp + as.numeric(object@phred[[i]])
            }
        }
    }
    mp <- matrix(mp, ncol = maxSeqLen)
    rownames(mp) <- rownames(object@phred[[1]])
    colnames(mp) <- 1:maxSeqLen
    return(mp)
})

setMethod("mergedPhredQuantiles", "Fastqq", function(object, quantiles)
{
    if(!(all(quantiles >= 0) & all(quantiles <= 1)) )
        stop("[mergedPhredQuantiles.Fastqq] all quantiles mustbe in [0, 1]")
    quantiles <- sort(unique(round(quantiles, 2)))
  
    sql <- seqLen(object)
    maxSeqLen <- max(sql[2, ])
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    ## Count qual values for each sequence position
    ## Convert counts into column-wise relative values
    ## Maximum counted read length
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    vec <- 1:maxSeqLen
    mrg <- mergedPhred(object)
    qrel <- as.data.frame(apply(mrg[, vec], 2, rel_real))
    names(qrel)
    names(qrel) <- vec
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    ## Walk through each column and extract row number
    ## for given quantile values
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    res <- .Call("get_column_quantiles", quantiles, qrel, PACKAGE="seqTools")
    return(res)
})


setMethod("plotMergedPhredQuant", "Fastqq", function(object, main, ...)
{
    quant <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    cols <- c("#1F78B4", "#FF7F00", "#E31A1C", "#FF7F00", "#1F78B4")
    qq <- mergedPhredQuantiles(object, quant)
    maxQ = floor(1.2*max(qq))
    xv <- 1:ncol(qq)
    
    if(missing(main))
        main <- paste("Merged position wise Phred Quantiles.", sep = "")
    
    plot(xv, xv, ylim=c(0, maxQ), type="n", bty="n", las=1,
        ylab="Phred score", xlab="Sequence position", main=main, ...)
    
    
    lines(xv, qq[1, ], col=cols[1], lty=2)
    lines(xv, qq[2, ], col=cols[2], lty=1)
    lines(xv, qq[3, ], col=cols[3], lwd=2)
    lines(xv, qq[4, ], col=cols[4], lty=1)
    lines(xv, qq[5, ], col=cols[5], lty=2)
    
    legend("top", ncol=6, lty=c(2, 1, 1, 1, 2), 
            lwd=c(1, 1, 2, 1, 1), col=cols, xjust=0.5, 
            legend=c("10%", "25%", "50%", "75%", "90%"), bty="n", cex=0.8)
    
    return(invisible()) 
})

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## End: Phred related functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Nucleotide frequency related functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##


setMethod("plotNucFreq", "Fastqq", function(object, i, main, maxx, ...)
{
    if(!is.numeric(i))
        stop("'i' must be numeric.")
    
    i <- as.integer(i)
    
    if( (i < 1) || (i > object@nFiles) )
        stop("'i' must be >0 and <nFiles.")
    
    maxSeqlen <- max(seqLen(object)[2, ])
    if(missing(maxx))
    {
        maxx <- maxSeqlen
    }
    else
    {
        if(!is.numeric(maxx))
            stop("'maxx' must be numeric.")
        maxx <- as.integer(maxx)
        if(maxx < 1)
            stop("'maxx' must be >0.")
        if(maxx > maxSeqlen)
            maxx <- maxSeqlen
    }
    
    # Skip extra line
    x <- 1:maxx
    nac <- object@nac[[i]][1:4, x]
    nacrel <- apply(nac, 2, rel_int)
    maxY = round(1.4 * max(nacrel), 1)
    
    # Maximum counted read length
    nacrel <- apply(nac, 2, rel_int)
    cols <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00")
    
    if(missing(main))
        main <- paste("Position wise Nucleotide frequency (",
                                probeLabel(object)[i], ")", sep="")
  
    plot(x, x, ylim=c(0, maxY), type="n", bty="n", las=1,
                ylab="Nucleotide fequency", xlab="sequence position",
                main=main, ...)
 
    lines(x, nacrel[1, ], col=cols[1], lwd=2)
    lines(x, nacrel[2, ], col=cols[2], lwd=2)
    lines(x, nacrel[3, ], col=cols[3], lwd=2)
    lines(x, nacrel[4, ], col=cols[4], lwd=2)
 
    legend("top", ncol=6, 
            lwd=c(1, 1, 2, 1, 1), col=cols, xjust=0.5, 
            legend=c("A", "C", "G", "T"), bty="n", cex=0.8)
    
    return(invisible())
})


setMethod("plotGCcontent", "Fastqq", function(object,main,...)
{
    cols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
            "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    ## Normalize matrix to colsum = 1
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    
    gc <- apply(object@gcContent, 2, rel_int)
    maxY = round(1.3*max(gc), 2)
    nf <- nFiles(object)
    x <- 1:nrow(gc)
    
    if(missing(main))
        main<-"GC content"
    
    plot(x, x, ylim=c(0, maxY), type="n", bty="n", las=1,
        ylab="Proportion of reads (%)", xlab="Relative GC content (%)", 
        main=main, ...) 
  
    for(i in 1:nf)
        lines(x, gc[, i], col=cols[i], lwd=2)
  
    legend("right", lwd=2, col=cols, legend=probeLabel(object),
                                                bty="n", cex=0.8)
  
    return(invisible())
})


setMethod("plotNucCount", "Fastqq", function(object, nucs=16, maxx, ...)
{
  
    # j = 16: N,    j = 2:3: gc
    if(!is.numeric(nucs))
        stop("'nucs' must be numeric.")
    
    nucs <- as.integer(nucs)
    if(any(nucs < 1) || any(nucs > 19))
        stop("'nucs' must be >0 and <20.")

    maxSeqlen <- max(seqLen(object)[2, ])
    
    if(missing(maxx))
    {
        maxx <- maxSeqlen
    }
    else
    {
        if(!is.numeric(maxx))
            stop("'maxx' must be numeric.")
        
        maxx <- as.integer(maxx)
        
        if(maxx<1)
            stop("'maxx must be >0.")
        
        if(maxx > maxSeqlen)
            maxx <- maxSeqlen
    }
    
    n <- nFiles(object)
    fvec <- 1:n
    
    i <- 1
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Skip extra line
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    x <- 1:maxx
    nac <- object@nac[[i]][, x]
    nacrel <- apply(nac, 2, rel_int)
    
    if(length(nucs) == 1){
        dfr <- data.frame(a = nacrel[nucs, ])
    }else{
        dfr <- data.frame(a = apply(nacrel[nucs, ], 2, sum))
    }
    
    if(n > 1)
    {
        for(i in 2:n)
        {
            nac <- object@nac[[i]][, x]
            nacrel <- apply(nac, 2, rel_int)
            if(length(nucs) == 1){
                dfr[, i] <- nacrel[nucs, ]
            }else{
                dfr[, i] <- apply(nacrel[nucs, ], 2, sum)
            }
        }        
    }
    
    maxY = round(1.4 * max(dfr), 3)
    cols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
    
    nv <- paste(iupac_chars[nucs], collapse = "")
    
    plot(x, x, ylim=c(0, maxY), type="n", bty="n", las=1,
        ylab="Nucleotide fequency", xlab="sequence position",
        main=paste("Position wise Nucleotide frequency:  '",
                                            nv, "'", sep=""))
  
    for(i in fvec)
        lines(x, dfr[, i], col=cols[i %% 10], lwd=2)
    
    legend("top", ncol=6, 
        lwd=c(1, 1, 2, 1, 1), col=cols[fvec %% 10], xjust=0.5,
        legend=probeLabel(object), bty="n", cex=0.8)
    
    return(invisible())
})

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## End: Nucleotide frequency related functions
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##



setMethod("plotKmerCount", "Fastqq",
                    function(object, index, mxey, main="K-mer count", ...)
{
    n <- nFiles(object)
    if(missing(index))
    {
        index <- 1:n
    }else{
        if(!is.numeric(index))
            stop("'index' must be numeric.")
        
        index <- sort(unique(as.integer(index)))
        
        if(any(index < 0) || any(index > n))
                stop("'index' must be in 1, .., ", n, " ( = nFiles).")
    }
    
    if(!missing(mxey))
    {
        if(!is.numeric(mxey))
            stop("'mxey' must be numeric.")
        mxey <- as.integer(mxey)
        
        if(mxey<1)
            stop("'mxey' must be positive.")
    }
    
    cols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")  
    pk <- 6
    if(object@k <= pk)
    {
        pk <- object@k
        x <- 1:(4^pk)
        
        if(missing(mxey))
            lg2y <- floor(log2(1.2 * (max(object@kmer))))
        else
            lg2y <- mxey
        
        maxY <- 2^lg2y

        
        plot(x, x, ylim=c(0, maxY), type="n", bty="n", las=1,
             ylab="K-mer count", xlab="K-mer index", main=main,
                                                axes=FALSE, ...)
    
    for(i in index)
        lines(x, object@kmer[, i], col=cols[i], lwd=2)
    }
    else
    {
        x <- 1:(4^pk)
        melt_factor <- as.integer(4^(object@k - pk))
        
        y_factor <- max(.Call("melt_vector", object@kmer[, 1], melt_factor,
                        PACKAGE="seqTools")) / max(object@kmer[, 1])
    
        if(missing(mxey))
            lg2y <- floor(log2(1.2 * (max(object@kmer)) * y_factor))
        else
            lg2y <- mxey
        maxY <- 2^lg2y
    
        plot(x, x, ylim=c(0, maxY), type="n", bty="n", las=1,
            ylab="K-mer count", xlab="K-mer index", 
                            main=main, axes=FALSE, ...)
    
        for(i in index)
        {
            lines(x, .Call("melt_vector", object@kmer[, i], melt_factor,
                        PACKAGE="seqTools"), col=cols[i], lwd=2)
        }
    }
    axis(side=1, at=0:4 * 4^(pk - 1), labels=c("A", "C", "G", "T", ""))
  
    axis(side=2, at=c(0, maxY),
                    labels=c(0, paste("2^", lg2y, sep="")), las=1)
  
    legend("right", lwd=2, col=cols[index], 
                legend=probeLabel(object)[index], bty="n", cex=0.8)
  
    return(invisible())
})



## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Merging and melting Fastqq objects
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
setMethod("mergeFastqq", "Fastqq", function(lhs, rhs)
{
    
    res <- new("Fastqq")
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    ## Simple concatenations
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    
    res@filenames <- c(lhs@filenames, rhs@filenames)
    res@nFiles <- lhs@nFiles+rhs@nFiles
    res@nReads <- c(lhs@nReads, rhs@nReads)
    res@nN <- c(lhs@nN, rhs@nN)
    res@seqLenCount <- cbind(lhs@seqLenCount, rhs@seqLenCount)
    res@gcContent <- cbind(gcContentMatrix(lhs), gcContentMatrix(rhs))
    res@seqLen <- cbind(lhs@seqLen, rhs@seqLen)
    
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    ## Singularize probeLabel
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    
    res@probeLabel <- c(lhs@probeLabel, rhs@probeLabel)  
    if(any(table(res@probeLabel) > 1))
    {
        message("[mergeFastqq] Singularizing probeLabel (new suffix).")
        res@probeLabel <- paste(1:res@nFiles, res@probeLabel, sep="_")
    }
    
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    ## Eventually resize arrays when lhs and rhs have different maxSeqLen
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
    if(lhs@maxSeqLen > rhs@maxSeqLen)
    {
        message("[mergeFastqq] Resizing rhs.")
        msl <- lhs@maxSeqLen
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        ## nac
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
        new_dim <- as.integer(c(nrow(rhs@nac), msl))
        rhs_nac <- .Call("r_enlarge_int_mat", rhs@nac, new_dim, 
                                                    PACKAGE="seqTools")
        res@nac <- c(lhs@nac, rhs_nac)
        
        # seqLenCount
        new_dim <- as.integer(c(msl, rhs@nFiles))
        rhs_seqLenCount <- .Call("r_enlarge_int_mat", 
                                rhs@seqLenCount, new_dim, PACKAGE="seqTools")
        
        res@seqLenCount <- c(lhs@seqLenCount, rhs_seqLenCount)
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        ## phred
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
        new_dim <- as.integer(nrow(rhs@phred), msl)
        rhs_phred_list <- list()
        
        for(i in 1:nFiles(rhs)){
                rhs_phred_list[[i]] <- .Call("r_enlarge_int_mat", 
                                rhs@phred[[i]], new_dim, PACKAGE="seqTools")
        }
        
        res@phred <- c(lhs@phred, rhs_phred_list)
        res@maxSeqLen <- lhs@maxSeqLen
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
    } else if(rhs@maxSeqLen > lhs@maxSeqLen)
    {
        
        message("[mergeFastqq] Resizing lhs.")
        msl <- rhs@maxSeqLen
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        ## nac
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
        new_dim <- as.integer(c(nrow(lhs@nac), msl))
        
        lhs_nac <- .Call("r_enlarge_int_mat",
                                lhs@nac, new_dim, PACKAGE="seqTools")
        
        res@nac <- c(rhs@nac, lhs_nac)
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        ## seqLenCount
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
        new_dim <- as.integer(c(msl, lhs@nFiles))
        
        lhs_seqLenCount <- .Call("r_enlarge_int_mat",
                    lhs@seqLenCount, new_dim, PACKAGE="seqTools")
        
        res@seqLenCount <- c(rhs@seqLenCount, lhs_seqLenCount)
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        ## phred
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
        new_dim <- as.integer(nrow(lhs@phred), msl)
        lhs_phred_list <- list()
        
        for(i in 1:nFiles(lhs))
        {
            lhs_phred_list[[i]] <- .Call("r_enlarge_int_mat", 
                    lhs@phred[[i]], new_dim, PACKAGE="seqTools")
        }
        res@phred <- c(rhs@phred, lhs_phred_list)
        res@maxSeqLen <- rhs@maxSeqLen
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
    } else { 
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        ## rhs@maxSeqLen == lhs@maxSeqLen
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
        res@maxSeqLen = lhs@maxSeqLen
        res@nac <- c(lhs@nac, rhs@nac)
        res@phred <- c(lhs@phred, rhs@phred)
        
        ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
        
    }
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Eventually melt down k-mer count matrix
    ## when lhs and rhs have different k
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    res@k <- pmin(lhs@k, rhs@k)
    
    kml <- kmerCount(lhs)
    if(lhs@k > res@k)
    {
        kml <- .Call("melt_kmer_matrix",
                    kml, c(lhs@k, res@k), PACKAGE="seqTools")
    }
    
    kmr <- kmerCount(rhs)
    if(rhs@k > res@k)
    {
        kmr <- .Call("melt_kmer_matrix",
                    kmr, c(rhs@k, res@k), PACKAGE="seqTools")
    }
    res@kmer <- cbind(kml, kmr)
    
    fkml <- lhs@firstKmer
    if(lhs@k > res@k)
    {
        fkml <- .Call("melt_kmer_matrix",
                    fkml, c(lhs@k, res@k), PACKAGE="seqTools")
    }
    fkmr <- rhs@firstKmer
    
    if(rhs@k > res@k)
    {
        fkmr <- .Call("melt_kmer_matrix",
                    fkmr, c(rhs@k, res@k), PACKAGE="seqTools")
    }
    res@firstKmer <- cbind(fkml, fkmr)
    
    return(res)
})


setMethod("meltDownK", "Fastqq", function(object, newK)
{
  
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    if(!is.numeric(newK))
        stop("'newK' must be numeric")
    
    newK <- as.integer(newK)
    
    if(length(newK) != 1)
        stop("'newK' must have length 1.")
    
    if(newK < 1 || newK >=  getK(object))
        stop("'newK' must be >= 1 and <= ", getK(object), ".")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    res <- new("Fastqq")
    res@filenames <- object@filenames
    res@nFiles <- object@nFiles
    res@k <- newK
    res@maxSeqLen <- object@maxSeqLen
    
    res@kmer <- .Call("melt_kmer_matrix",
        object@kmer, c(getK(object), newK), PACKAGE="seqTools")
    
    res@firstKmer <- .Call("melt_kmer_matrix",
        object@firstKmer, c(getK(object), newK), PACKAGE="seqTools")
    
    res@nReads <- object@nReads
    res@nN <- object@nN
    res@seqLenCount <- object@seqLenCount
    res@gcContent <- object@gcContent
    res@nac <- object@nac
    res@phred <- object@phred
    res@seqLen <- object@seqLen
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
  
    return(res)
})


listMelt <- function(x, oldK, newK)
{
    f <- function(x) .Call("melt_kmer_matrix", x, 
            as.integer(c(oldK, newK)), PACKAGE="seqTools")
  
    return(lapply(x, f))
}



setMethod("[", "Fastqq", function(x, i, j, drop="missing")
{
  
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    res <- new("Fastqq")
    res@filenames <- x@filenames[i]
    res@probeLabel <- x@probeLabel[i]
    res@nFiles <- length(i)
    res@k <- x@k
    res@seqLen <- x@seqLen[, i, drop=FALSE]
    res@maxSeqLen <- max(res@seqLen[2, ])
    res@kmer <- x@kmer[, i, drop=FALSE]
    res@firstKmer <- x@firstKmer[, i, drop=FALSE]
    res@nN <- x@nN[i]
    res@seqLenCount <- x@seqLenCount[, i, drop=FALSE]
    res@gcContent <- x@gcContent[, i, drop=FALSE]
    res@nac <- x@nac[i]
    res@phred <- x@phred[i]
    res@nReads <- x@nReads[i]
    res@collectTime <- x@collectTime
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    
    return(res)
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Calculation of distance matrix based on Canberra distance
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

setMethod("cbDistMatrix", "Fastqq", 
        function(object, nReadNorm=max(nReads(object)))
{
    if(!is.numeric(nReadNorm))
        stop("'nReadNorm' must be numeric.")
    nReadNorm <- as.integer(nReadNorm)
    
    if(nReadNorm < max(nReads(object)))
        stop("'nReadNorm' must be greater than all nRead.")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Column-wise normalizing read counts (by upscaling) 
    ## so that column sums become nearly equal in order to
    ## compensate sequencing depth artifacts in 
    ## Canberra distance values
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    scale <- nReadNorm/nReads(object)
    
    scaled <- .Call("scale_kmer_matrix",
            kmerCount(object), scale, PACKAGE="seqTools")
    
    colnames(scaled) <- object@probeLabel
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
  
    return(.Call("cb_distance_matrix", scaled, PACKAGE="seqTools"))
})


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
