/*
 * fq_traverse.h
 *
 *  Created on: 12.10.2013
 *      Author: wolfgang
 */

#ifndef FQ_TRAVERSE_H_
#define FQ_TRAVERSE_H_

#include "stat_defs.h"
#include "dna_astream.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 	Fastq traverse:
//				Contains dna_astream which is traversed.
//				The <seqname> following '+' is optional,
//					but if it appears right after '+', it should be identical to the <seqname> following '@'.
//				The length of <seq> is identical the length of <qual>. Each character in <qual> represents the phred quality of the corresponding nucleotide in <seq>.

//
///////////////////////////////////////////////////////////////////////////////////////////////////

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Definition of static constants
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static const unsigned long  faq_array_size  	=0xFFFFFF;	// char array capacity

static const char			faq_char_sqDelim 	='@';	// Fastq sequence delimiter
static const char			faq_char_qlDelim 	='+';	// Fastq quality  delimiter
static const char			faq_char_lf			='\n';	// Line feed
static const char			faq_char_eos		='\0';	// End of string

static const int			faq_ok				=0;
static const int			faq_empty			=1;
static const int			faq_newSeq			=4;


typedef struct fastq_traverse
{
	daStream *das;

	int state;
	unsigned nSeq;		// Number of sequences
	unsigned nFill;		// Number of fill operations.
	unsigned nProcFull;	// Number of times where proc array was filled up.

	unsigned minSeqLen;
	unsigned maxSeqLen;

	unsigned lastSeqLen;
	unsigned nUneqLeqLen;

} fqTraverse;


int faqEmpty(fqTraverse* faq)	{ return dasEmpty(faq->das);  }
int faqIsOpen(fqTraverse* faq)	{ return dasIsOpen(faq->das); }
int faqIsEof(fqTraverse* faq)	{ return dasIsEof(faq->das); }

void faq_destroy(fqTraverse *faq)
{
	if(faq)
	{
		das_destroy(faq->das);
		free(faq);
	}
}


fqTraverse * faq_init(const char* filename, unsigned mode)
{
	fqTraverse * faq=(fqTraverse*)calloc(sizeof(fqTraverse),1);

	if(!faq)
		return 0;

	faq->das=das_init(filename,mode,faq_array_size);

	if(!faq->das)
	{
		//printf("[faq_init] das_init returned 0!\n");
		free(faq);
		return 0;
	}
	// initialize with large value
	--faq->minSeqLen;

	// Returns memory initialized object.
	// File may be closed.
	// File may be empty.
	// Stream is not yet initialized.
	return faq;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Check-Fill: Checks for empty rfc (frs==0) and eventually fills das.
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

int faqCheckFill(fqTraverse *faq)
{
	//printf("[faqCheckFill].\n");
	if(dasEmpty(faq->das))
	{
		//printf("[faqCheckFill] ->das_fill\n");
		if(das_fill(faq->das))
		{
			++(faq->nFill);
			faq->state|=faq_empty;
			return faq->state;
		}
		faq->state&=(~faq_empty);
	}
	return faq_ok;
}

int faqCheckCapFill(fqTraverse *faq,int capacity)
{
	daStream *das=faq->das;
	if(das->r_end-das->r_iter<capacity)
	{
		if(das_fill(das))
		{
			++(faq->nFill);
			faq->state|=faq_empty;
			return faq->state;
		}
		faq->state&=(~faq_empty);
	}
	return faq_ok;
}



// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Check section (Must be placed beforehand skip-routines)
//		faq_checkNewLine
//		faq_checkSeqHeader
//		faq_checkComment
//		faq_checkN
//		faq_check_nuc
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //



// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Skip section (Do not reorder...)
//			faq_checkSkipNewLine
//			faq_skipLine
//			faq_skipComment
//			faq_skipN
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static inline int faq_checkSkipNewLine(fqTraverse *faq)
{
	// Skips LF ('\n') character and
	// eventually tries to fill array

	//printf("[CheckNewLine].\n");
	if(*faq->das->r_iter==faq_char_lf)
	{
		if(faqCheckCapFill(faq,2))
			return faq->state;
		// Eat newline
		++faq->das->r_iter;
		return faq_ok;
	}
	return faq->state;
}

static inline int faq_skipLastNewLine(fqTraverse *faq)
{
	daStream *das=faq->das;
	if((*faq->das->r_iter==faq_char_lf) & (das->dnaf->state==dfs_file_closed))
	{
		++das->r_iter;
		return faq_ok;
	}
	return faq_empty;
}


static inline int faq_skipLine(fqTraverse *faq)
{
	// Traverses until first position of next line
	// and eventually tries to refill array
	daStream *das=faq->das;
	while(*das->r_iter!=faq_char_lf)
	{
		++das->r_iter;
		if(faqCheckFill(faq))
			return faq->state;
	}
	++das->r_iter;

	if(faqCheckFill(faq))
		return faq->state;

	return faq_ok;
}

int faq_procNuc(fqTraverse *faq)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Copies one line of continuous Nucleotide sequence
	// from:	input  fastq stream			(rfc array)
	// to  :	output nucleotide stream	(pos array)
	//
	// Copy procedure ends when:
	//	A) A non nucleotide is found.
	//	B) Quality header ('+') is found.
	//	C) rfc array runs empty.
	//	D) pos array is filled up.
	//
	// Tries to recover when end of rfc-array is reached.
	// Copying always starts at begin of nuc array
	//
	// Returns fat_ok (=0) when *NO* output is created.
	//
	// Allows 'if(fat_procNuc) => process output'
	// testing.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	daStream *das=faq->das;
	unsigned seqLen=0;

	//printf("[faq pNuc] start: '%s'\n",das->r_iter);

	// Empty array must be prevented here
	if(faqCheckFill(faq))
	{
		faq->lastSeqLen=seqLen;
		return seqLen;
	}


	if(*faq->das->r_iter==faq_char_sqDelim)
	{
		// Skip header line
		if(faq_skipLine(faq))
		{
			faq->lastSeqLen=seqLen;
			return seqLen;
		}

		// Goto begin of pos array.
		das->p_iter=das->pos;
		//printf("[faq pNuc] while: '%s'\n",das->r_iter);
		++faq->nSeq;

		while(!((ACGT[(unsigned)*das->r_iter]==zval) || *das->r_iter==faq_char_qlDelim ) )
		{

			*das->p_iter=*das->r_iter;
			++das->r_iter;
			++das->p_iter;
			++seqLen;
			//printf("[faq pNuc]  fill: '%s'\n",das->r_iter);
			// Two recoverable conditions:
			if(faqCheckFill(faq))
			{
				faq->lastSeqLen=seqLen;
				return seqLen;
			}
			if(faq_checkSkipNewLine(faq))
			{
				faq->lastSeqLen=seqLen;
				return seqLen;
			}
			if(dasProcEmpty(das))
			{
				++faq->nProcFull;
				faq->lastSeqLen=seqLen;
				return seqLen;
			}
		}
		// Add string termination
		if(!dasProcEmpty(das))
		{
			das->npPos=(das->p_iter-das->pos);
			*das->p_iter=faq_char_eos;
		}
	}
	faq->lastSeqLen=seqLen;
	return seqLen;
}

int faq_procPhred(fqTraverse *faq)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Copies one line of Phred characters
	// from:	input  fastq stream			(rfc array)
	// to  :	output nucleotide stream	(pos array)
	//
	// Copy procedure ends when:
	//	A) Sequence header ('@') is found.
	//	B) rfc  array runs empty.
	//	C) pos array is filled up.
	//
	// Tries to recover when end of rfc-array is reached.
	// Copying always starts at begin of pos array.
	//
	// Returns fat_ok (=0) when *NO* output is created.
	//
	// Allows 'if(fat_procNuc) => process output'
	// testing.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	daStream *das=faq->das;
	unsigned seqLen=0;

	// Empty array must be prevented here
	if(faqCheckFill(faq))
	{
		faq->lastSeqLen=seqLen;
		return seqLen;
	}

	//printf("[faq_Phred] start: '%s'\n",das->r_iter);
	if(*faq->das->r_iter==faq_char_qlDelim)
	{
		// Skip header line
		if(faq_skipLine(faq))
		{
			faq->lastSeqLen=seqLen;
			return seqLen;
		}
		if(faqCheckFill(faq))
			return 0;

		// Goto begin of pos array.
		das->p_iter=das->pos;
		// Copy until new seq header is found

		//printf("[faq_Phred] while: '%s'\n",das->r_iter);
		while(!(*das->r_iter==faq_char_sqDelim))
		{
			*das->p_iter=*das->r_iter;
			++das->r_iter;
			++das->p_iter;
			++seqLen;

			//printf("[faq Phred]  while: '%s'\tr_iter:%li\n",das->r_iter,das->r_end-das->r_iter);
			// Two recoverable conditions:
			if(faqCheckFill(faq))
			{
				//printf("[faq_Phred] checkFill: %u\n",faq->nFill);
				if(seqLen!=faq->lastSeqLen)
					++faq->nUneqLeqLen;
				return seqLen;
			}
			if(faq_checkSkipNewLine(faq))
			{
				if(seqLen!=faq->lastSeqLen)
					++faq->nUneqLeqLen;
				return seqLen;
			}
			if(dasProcEmpty(das))
			{
				++faq->nProcFull;
				if(seqLen!=faq->lastSeqLen)
					++faq->nUneqLeqLen;
				return seqLen;
			}
		}
		// Add string termination
		if(!dasProcEmpty(das))
		{
			das->npPos=(das->p_iter-das->pos);
			*das->p_iter=faq_char_eos;
		}

		if(seqLen!=faq->lastSeqLen)
			++faq->nUneqLeqLen;
	}
	return seqLen;
}


#endif /* FQ_TRAVERSE_H_ */
