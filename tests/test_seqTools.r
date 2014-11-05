
library(seqTools)

filename <- "test_seqTools.r"
basedir <- system.file("extdata", package = "seqTools")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## countDnaKmers
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

if(!identical(countDnaKmers("ACGT", k = 1, start = 3:1, width = 1), cdk_ACGT))
    stop("[countDnaKmers] Test 1 '", filename, "' FAILED!")
   
if(!identical(countDnaKmers("ACGT", k = 1, start = 3, width = 1), cdk_ACGT_one))
    stop("[countDnaKmers] Test 2 '", filename, "' FAILED!")

if(!identical(
        countDnaKmers("ATTNAC", k = 2, start = 1:3, width = 1), cdk_ATTNAC))
    stop("[countDnaKmers] Test 3 '", filename, "' FAILED!")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## revCountDnaKmers
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

if(!identical(
        revCountDnaKmers("ACGTACGT", k = 2, start = 6:4, width = 2), rck_ACGT))
    stop("[revCountDnaKmers] Test 1 '", filename, "' FAILED!")


## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## kmerCount.fastqq
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

fq <- fastqq(file.path(basedir, "test_l5_N.fq"), k = 2)
if(!identical(kmerCount(fq), kmer_l5_N))
    stop("[kmerCount.fastqq] Test 1 '", filename, "' FAILED!")

fq<-fastqq(file.path(basedir, "test_l6.fq"), k = 2)
if(!identical(kmerCount(fq), kmer_l6))
    stop("[kmerCount.fastqq] Test 2 '", filename, "' FAILED!")

fq<-fastqq(file.path(basedir, "test_l6_multi_line.fq"), k = 2)
if(!identical(kmerCount(fq), kmer_l6_ml))
    stop("[kmerCount.fastqq] Test 3 '", filename, "' FAILED!")

fq<-fastqq(file.path(basedir, "test_l10_20_40.fq"),k = 2)
if(!identical(kmerCount(fq), kmer_l10_20))
    stop("[kmerCount.fastqq] Test 4 '", filename, "' FAILED!")

fq<-fastqq(file.path(basedir, "test_l10_atcg.fq"), k = 2)
if(!identical(kmerCount(fq), kmer_l10_atcg))
    stop("[kmerCount.fastqq] Test 5 '", filename, "' FAILED!")

fq<-fastqq(file.path(basedir, "test_l10_ATCGN.fq"), k = 2)
if(!identical(kmerCount(fq), kmer_l10_ATCGN))
    stop("[kmerCount.fastqq] Test 6 '", filename, "' FAILED!")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## ascii2char, char2ascii
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

if(!identical(ascii2char(97:101, multiple = FALSE), "abcde"))
    stop("[ascii2char] Test 1 '", filename, "' FAILED!")

if(!identical(ascii2char(97:101, multiple = TRUE), letters[1:5]))
    stop("[ascii2char] Test 2 '", filename, "' FAILED!")

if(!identical(ascii2char(char2ascii("abcde")), "abcde"))
    stop("[ascii2char] Test 3 '", filename, "' FAILED!")

if(!identical(char2ascii("abcde"), 97:101))
    stop("[char2ascii] Test 1 '", filename, "' FAILED!")

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## END OF FILE
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
