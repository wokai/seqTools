/*
 * seqTools.h
 *
 *  Created on: 14.10.2013
 *      Author: kaisers
 *
 *	C99 standard: GCC C compiler / Debugging / Other debugging flags: -std=c99
 *
 * Project	:	seqTools
 * Created	:	26.08.2013
 * Author	:	W. Kaisers
 *
 * Content	:	Counting and evaluation of DNA-Motifs (short DNA sequences of length 5-10)
 *
 * Version	:	0.99.7
 *
 * Changelog	:
 * 14.10.2013	:	Creation of Project
 * 04.12.2013	:	Submission to Bioconductor
 * 15.07.2014	:	Successfully linked to zlibbioc
 * 24.07.2014   :       Added news file
 *
 */

#ifndef SEQTOOLS_H_
#define SEQTOOLS_H_

// /usr/share/R/include

#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <locale.h>			// thousands separator
#include <stdbool.h>		// bool (true, false)
#include <math.h>  			// log2 (-lm)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h> // DllInfo
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>

#include "stat_defs.h"
#include "fa_traverse.h"
#include "fastq_parser.h"
#include "resize_matrix.h"
#include "seq_art.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		DNA k-mers
//
///////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * 		k			:	number of nucleotides in DNA motif (Kmer length)
 * 						4^n possible motifs for n = mset
 *
 * 		mindex		:	0-based unique index for each motif in mset.
 * 						Takes values in { 0, ..., (4^n-1) }
 *
 * + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
 *
 *	Window - Frame nomenclature:
 *	Window : Range which contains all (1-based) positions for Kmer-starts
 *	Frame  : Range which contains all (1-based) positions which belong to some nMier in window
 *			 Frame overhangs window for (n-1)   positions in reading downstream direction
 *
 * 						               10
 * 						1     5       9      3  6
 * 						TTTT | CCCC | GGGG | AAAA
 * 	 n=3
 * 	 + strand                  wwww---ff
 * 	 x =start                  x-->
 *
 *   - strand                           ff---wwww
 *   x =start                                <--x
 *
 * + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
 *
 *
 */



///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		count_Kmers
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP count_Kmers(SEXP pSeq, SEXP pK, SEXP pWidth, SEXP pNn);
SEXP count_dna_Kmers(SEXP pSeq, SEXP pStart, SEXP pWidth, SEXP pK, SEXP pNn);
SEXP count_genome_Kmers(SEXP pSeq, SEXP pSeqid, SEXP pLstart, SEXP pWidth, SEXP pStrand, SEXP pK, SEXP pNn);
SEXP count_splice_Kmers(SEXP pSeq, SEXP pSeqid, SEXP pLend, SEXP pRstart, SEXP pWidth, SEXP pStrand, SEXP pK, SEXP pNn);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		count_fastq
//
///////////////////////////////////////////////////////////////////////////////////////////////////

const int default_fastq_max_seqlen = 200;
const int default_fq_buf_capacity = 10000;

SEXP count_fastq(SEXP pInfile, SEXP pK);
SEXP trim_fastq(SEXP pInfile, SEXP pVals, SEXP pOutfile);
SEXP fastq_KmerSubset_locs(SEXP pInfile, SEXP pK, SEXP pKmerIdx);
SEXP get_Kmer_Index(SEXP pSequence, SEXP pK);
SEXP get_kmer(SEXP pKmerIndex, SEXP pK);
SEXP fastq_Kmer_locs(SEXP pInfile, SEXP pK);

const unsigned default_fasta_Kmers_col_number=50;

SEXP write_fai(SEXP pInfile, SEXP pOutfile);
SEXP count_fasta_Kmers(SEXP pFasta, SEXP pK);


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		rev_count_Kmers
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP rev_count_dna_Kmers(SEXP pSeq, SEXP pStart, SEXP pWidth, SEXP pK, SEXP pNn);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		get_column_quantiles
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// Expects a vector of quantiles (pQuant) and
// data.frame (pDf) where each column contains relative quantities (sums up to 1)
// Steps down each column and writes values for each quantile
// into output data.frame

SEXP get_column_quantiles(SEXP pQuant, SEXP pDf);
SEXP melt_vector(SEXP pValues, SEXP pFactor);
SEXP melt_kmer_matrix(SEXP pKmerCount, SEXP pK);
SEXP scale_kmer_matrix(SEXP pKmerCount, SEXP pScale);
SEXP cb_distance_matrix(SEXP pKmerCount);
SEXP rel_int(SEXP pInt);
SEXP rel_real(SEXP pReal);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simulation of fastq Reads
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP sim_k_values(SEXP pVal);
SEXP sim_dna_k_mer(SEXP pVal);
SEXP set_dna_k_mer(SEXP pVal, SEXP pSeq, SEXP pPosition, SEXP pKmerIndex, SEXP pSeqIndex);
SEXP gzwrite_fastq_dna(SEXP pVal, SEXP pSeq, SEXP pFilename);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations for R_registerRoutines
//
///////////////////////////////////////////////////////////////////////////////////////////////////

void R_init_seqTools(DllInfo *info);



#endif /* SEQTOOLS_H_ */
