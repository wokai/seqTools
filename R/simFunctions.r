
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
##                                                                           ##
##  Project   :   seqTools                                                   ##
##  Created   :   26.August.2013                                             ##
##  Author    :   W. Kaisers                                                 ##
##  Version   :   0.99.13                                                    ##
##  File      :   simFunctions.r                                             ##
##  Content   :   Functions which produce fastq files with simulated data    ##
##                (Done for testing Fastqq and DNA k-mer counting)           ##
##                writeSimFastq, writeSimContFastq, simFastqqRunTimes,       ##
##                sim_fq                                                     ##
##                                                                           ##
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Create fastq files with simulated k-mers
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

writeSimFastq<-function(k=6, nk=5, nSeq=10, filename="sim.fq.gz")
{
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # k
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(k))
        stop("'k' must be numeric.")
    if(length(k) != 1)
        stop("'k' must have length 1.")
    k <- as.integer(k)    
    if( (k < 1) || (k > max_k) )
        stop("'k' out of range.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # nk
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(nk))
        stop("'nk' must be numeric.")
    if(length(nk) != 1)
        stop("'nk' must have length 1.")
    nk<-as.integer(nk)
    if( (nk < 1) || (nk > 1000) )
        stop("nk must be positive and <1000.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # nSeq
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(nSeq))
        stop("'nSeq' must be numeric.")
    if(nSeq < 1)
        stop("'nSeq' must be positive.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # filename
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.character(filename))
        stop("'filename' must be string.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Do the work
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    bm <- Sys.localeconv()[7]
    val <- as.integer(c(k, nk, nSeq))

    pseq <- .Call("sim_dna_k_mer", val, PACKAGE="seqTools")
    res <- .Call("gzwrite_fastq_dna", val, pseq, filename, PACKAGE="seqTools")

    message("[writeSimFastq] file '", basename(filename), "': ",
            format(res,big.mark=bm), " Bytes written.")
  
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Terminate function.
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    return(invisible())
}

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## writeSimContFastq
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
writeSimContFastq<-function(k=6, nk=5, nSeq=10, pos=1, 
                kIndex=1, nContam=nSeq, filename="simc.fq.gz")
{
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # k
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(k))
        stop("k must be numeric.")
    if(length(k) != 1)
        stop("k must have length 1.")
    k<-as.integer(k)    
    if( (k < 1) || (k > max_k) )
        stop("k out of range.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # nk
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(nk))
        stop("nk must be numeric.")
    if(length(nk) != 1)
        stop("nk must have length 1.")
    nk<-as.integer(nk)
    if( (nk < 1) || (nk > 1000) )
        stop("nk must be positive and <1000.") 
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # nSeq
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(nSeq))
        stop("nSeq must be numeric.")
    if(length(nSeq) != 1)
        stop("nSeq must have length 1.")
    if(nSeq < 1)
        stop("nSeq must be positive.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # pos
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(pos))
        stop("pos must be numeric.")
    pos<-as.integer(pos)
    if( any(pos < 1) || any(pos > (k * (nk - 1) + 1)))
        stop("pos out of range.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # writeSimContFastq should take 1-based positions
    # while "set_dna_k_mer" takes    0-based positions
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    pos<-pos - 1L
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # kIndex
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(kIndex))
        stop("kIndex must be numeric.")
    
    kIndex <- as.integer(kIndex)
    
    if(length(pos) != length(kIndex))
        stop("pos and kIndex must have equal length.")
    
    if( any(kIndex < 1L) || any(kIndex > 4^k) )
        stop("kIndex out of range.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # writeSimContFastq should take 1-based k-indices
    # while "set_dna_k_mer" takes    0-based k-indices
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    kIndex <- kIndex - 1L
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # nContam: Deterministic replacements are copied
    # into the first 'nContam' reads
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.numeric(nContam))
        stop("nContam must be integer.")
    
    nContam <- as.integer(nContam)
    
    if(length(nContam) != 1)
        stop("nContam must have length 1.")
    
    if( (nContam < 1) || (nContam > nSeq) )
        stop("nContam must be >0 and <nSeq.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # filename
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    if(!is.character(filename))
        stop("filename must be character.")
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Do the work
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    bm <- Sys.localeconv()[7]
    
    val <- as.integer(c(k, nk, nSeq))
    
    pseq <- .Call("sim_dna_k_mer", val, PACKAGE="seqTools")
    
    prseq <- .Call("set_dna_k_mer", 
                val, pseq,pos, kIndex, nContam, PACKAGE="seqTools")
  
    res <- .Call("gzwrite_fastq_dna",
                    val, prseq, filename, PACKAGE="seqTools")
  
    message("[writeSimContFastq] file '", basename(filename), 
            "': ", format(res,big.mark=bm), " Bytes written.")
  
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Terminate function
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    return(invisible())
}



simFastqqRunTimes<-function(k, nSeq, filedir=".")
{
    if(missing(k))
        k <- 2:(max_k)
    else
    {
        if(!is.numeric(k))
            stop("'k' must be numeric.")
        
        k <- sort(unique(as.integer(k)))

        if(any(k < 1) || any(k > max_k))
            stop("'k' out of range.")
    }
    
    if(missing(nSeq))
        nSeq <- as.integer(c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7))
    else {
        if(!is.numeric(nSeq))
            stop("'nSeq' must be numeric.")
        nSeq <- as.integer(nSeq)
        if(any(nSeq < 1) )
            stop("'nSeq' must be positive.")
    }
    
    if(!is.character(filedir))
        stop("'filedir' must be character.")
    
    if(!file.exists(filedir))
    {
        if(!dir.create(filedir))
            stop("Cannot create filedir '", filedir, "'.")
    }
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Write simulated fastq files
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    message("[simFastqqRunTimes] Creating fastq files:")
    fqFiles <- character(length(nSeq))
    bm <- Sys.localeconv()[7]
    
    for(j in 1:length(nSeq))
    {
        fqFiles[j] <- file.path(filedir,
            paste("sfqrt_nSeq", nSeq[j], ".fq.gz", sep=""))
        
        message("[simFastqqRunTimes] (",format(j , w=3),"/", length(nSeq),
            ") nSeq=", format(nSeq[j], big.mark=bm, w=12))
        
        # Allways equal sized reads: k=6, read length=102
        writeSimFastq(k=6, nk=17, nSeq=nSeq[j], filename=fqFiles[j])    
    }
    
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    # Prepare result table
    # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
    nSim <- length(k) * length(nSeq)
    m <- 1
    res<-data.frame(id=1:nSim, k=0, nSeq=0, runtime=0)
    
    for(i in 1:length(k))
    {
        for(j in 1:length(nSeq))
        {
            message("[simFastqqRunTimes] (",
                        format(m, width=3),"/", nSim, ") Fastqq run.")
            
            fq <- fastqq(fqFiles[j], k=k[i])
            res$k[m] <- k[i]
            res$nSeq[m] <- nSeq[j]
            res$runtime[m] <- collectDur(fq)
            m <- m + 1
        }
    }
    return(res)
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Encapsulates whole simulation studies
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

sim_fq<-function(nRep=2, nContamVec=c(100, 1000), grSize=20, nSeq=1e4,
                k=6, kIndex=1365, pos=20)
{
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Arguments
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    # nContam   : Absolute number of contaminated sequences
    # grSize    : Nr of fastq files in control and contamination group
    # nSeq      : Nr of reads per fastq file
    # k         : k-mer length used in fastqq function
    # kIndex    : k-mer index of contaminating sequence default is "CCCCCC"
    # pos       : Contamination is inserted at position in read sequence, 
    #               default is 20
  
    if(length(kIndex) != length(pos))
        stop("'kIndex' and 'pos' must have equal length.")
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Fixed internal values
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    # ksim        : k-mer sized used for simulation
    ksim <- 6
    # nksim     : number of k-mers per read (determines read length)
    nksim <- 17
    
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Write simulated fastq files for control group
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    # ctrl_vec=leaf labels for control group
    ctrl_vec <- 1:grSize
    ctrl_files <- paste("sfq_i", ctrl_vec, "_ctrl.fq.gz", sep = "")
    
    for(i in ctrl_vec){
        writeSimFastq(k=ksim, nk=nksim, nSeq=nSeq,
                                        filename = ctrl_files[i])
    }
    
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    ## Prepare values for contamination and result
    ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
    
    # cont_vec=leaf labels for contamination group
    # (must be nonoverlapping with ctrl_vec)
    cont_vec <- ctrl_vec + 100
    contamVecLen <- length(nContamVec)
    nSim = nRep * contamVecLen
    # Result data.frame
    res <- data.frame(id=1:nSim, nContam=0, rep=0, sum=NA)
    
    m <- 1
    for(i in 1:length(nContamVec))
    {
        nContam<-nContamVec[i]
        for(j in 1:nRep)
        {
            message("[sim_fq] Simulation (",
                                    format(m, width=3), "/", nSim, ").")
            
            cont_files <- paste("sfq_i", cont_vec, "_cont.fq.gz", sep="")
            
            for(i in 1:length(cont_vec))
                writeSimContFastq(k=ksim, nk=nksim, nSeq=nSeq, pos=pos,
                    nContam=nContam, kIndex=kIndex,
                    filename=cont_files[i])

            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            ## read Fastqq
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            ctrl_fqq <- fastqq(ctrl_files, k=k)
            probeLabel(ctrl_fqq) <- ctrl_vec
            cont_fqq <- fastqq(cont_files, k=k)
            probeLabel(cont_fqq) <- cont_vec
            
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            ## merge, dist, hclust
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            mrg <- mergeFastqq(ctrl_fqq, cont_fqq)
            mtx <- cbDistMatrix(mrg)
            hc <- hclust(as.dist(mtx))
            
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            ## Calculate results and write into res
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            res$nContam[m] <- nContam
            res$rep[m] <- j
            
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            ## Calculate how many labels of the contamination group
            ## are in the first half of leaf labels
            ## + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
            res$sum[m] <- sum(hc$order[ctrl_vec] > grSize)
            m <- m + 1    
        }
    }
    res$lat <- pmin(res$sum, grSize - res$sum) / grSize * 100
    res$perc <- res$nContam / nSeq * 100
    return(res)
}


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
