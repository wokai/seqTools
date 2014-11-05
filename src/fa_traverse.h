/*
 * fa_traverse.h
 *
 *  Created on: 08.10.2013
 *      Author: kaisers
 */

#ifndef FA_TRAVERSE_H_
#define FA_TRAVERSE_H_

#include "stat_defs.h"
#include "dna_astream.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 	Fasta traverse:
//				Contains dna_astream which is traversed.
//				Identified Nucleotide sequences of length >= k will be copied to
//				das.pos array.
//
//				Defines subroutines for skipping header-lines, N-array and comment lines.
//
///////////////////////////////////////////////////////////////////////////////////////////////////


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// 	DNA file stream:
//				Struct layer for reading text into char * buffer
//				from either uncompressed or compressed files.
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Definition of static constants
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static const unsigned  		fat_array_size  =0x1000;	// char array capacity

static const char			fat_char_newSeq	='>';		// Fasta sequence delimiter
static const char			fat_char_cmt	=';'; 		// Fasta comment  delimiter
static const char			fat_char_lf		='\n';		// Line feed
static const char			fat_char_eos	='\0';		// End of string

static const int			fat_ok			=0;
static const int			fat_empty		=1;
static const int			fat_loc			=2;
static const int			fat_newSeq		=4;
static const int			fat_nucReady	=8;

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// 	Struct definition, basic accessors; constructor and destructors:
//		struct faTraverse
//		fatEmpty, fatNewSeq, fatNucReady
//		fat_destroy, fat_init
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

typedef struct fasta_traverse
{
	daStream *das;

	int state;
	unsigned k;

	unsigned nN;		// Number of N's so far seen
	unsigned nSeq;		// Number of sequences
	unsigned nCom;		// Number of comment lines
	unsigned nFill;		// Number of fill operations.

} faTraverse;

int fatEmpty (faTraverse* fat)  { return dasEmpty (fat->das); }
int fatIsOpen(faTraverse* fat)	{ return dasIsOpen(fat->das); }
int fatIsEof (faTraverse* fat)	{ return dasIsEof (fat->das); }
int fatNewSeq(faTraverse *fat)  { return fat->state & fat_newSeq; }

// Only purposed for use in top-level:
// resets nucReady when called.
int fatNucReady(faTraverse *fat)
{
	if(fat->state & fat_nucReady)
	{
		fat->state &= (~fat_nucReady);
		return 1;
	}
	return 0;
}


void fat_destroy(faTraverse *fat)
{
	if(fat)
	{
		das_destroy(fat->das);
		free(fat);
	}
}

faTraverse * fat_init(const char* filename, int k)
{
	faTraverse * fat = (faTraverse*) calloc(sizeof(faTraverse), 1);
	if(!fat)
		return 0;

	fat->das = das_init(filename, fat_array_size);
	if(!fat->das)
	{
		//printf("[fat_init] das_init returned 0!\n");
		free(fat);
		return 0;
	}
	fat->k = (unsigned) k;

	// Returns memory initialized object.
	// File may be closed.
	// File may be empty.
	// Stream is not yet initialized.
	return fat;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Check-Fill: Checks for empty rfc (frs==0) and eventually fills das.
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

int fatCheckFill(faTraverse *fat)
{
	//if(fat->das->frs==0)
	if(dasEmpty(fat->das) )
	{
		//printf("[fatCheckFill] ...\n");
		if(das_fill(fat->das))
		{
			++(fat->nFill);
			fat->state |= fat_empty;
			//printf("[fatCheckFill] return empty.\n");
			return fat->state;
		}
		fat->state &= (~fat_empty);
	}
	return fat_ok;
}


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Check section (Must be placed beforehand skip-routines)
//		fat_checkNewLine
//		fat_checkNewSeq
//		fat_checkComment
//		fat_checkN
//		fat_check_nuc
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE int fat_checkNewLine(faTraverse *fat)
{
	if(LFCR[(unsigned)*fat->das->r_iter] == lfcv)
	{
		fat->state |= fat_loc;
		return fat->state;
	}
	return fat_ok;
}

static R_INLINE int fat_checkNewSeq(faTraverse *fat)
{
	if(*fat->das->r_iter == fat_char_newSeq)
	{
		fat->state |= fat_newSeq;
		return fat->state;
	}
	return fat_ok;
}

static R_INLINE int fat_checkComment(faTraverse *fat)
{
	if(*fat->das->r_iter==fat_char_cmt)
	{
		//printf("[checkCmt] Comment found!\n");
		fat->state|=fat_loc;
		return fat->state;
	}
	return fat_ok;
}

static R_INLINE int fat_checkN(faTraverse *fat)
{
	if(ACGTN[(unsigned)*fat->das->r_iter] == nval)
	{

		fat->state |= fat_loc;
		return fat->state;
	}
	return fat_ok;
}

static R_INLINE int fat_check_nuc(faTraverse *fat)
{
	if(ACGT[ ((unsigned) *fat->das->r_iter) ] == zval)
	{
		fat->state |= fat_loc;
		return fat->state;
	}
	return fat_ok;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Skip section (Do not reorder...)
//			fat_skipNewLine
//			fat_skipLine
//			fat_skipComment
//			fat_skipN
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE int fat_skipNewLine(faTraverse *fat)
{
	// Skips LF ('\n', ASCII 10)  and CR (ASCII 13) character and
	// eventually tries to fill array
	if(LFCR[ (unsigned)*fat->das->r_iter ] == lfcv)
	{
		if(fatCheckFill(fat))
			return fat->state;

		// Eat newline
		++fat->das->r_iter;

		// Unset flag
		fat->state &= (~fat_loc);

		return fat_ok;
	}
	return fat->state;
}

static R_INLINE int fat_skipLine(faTraverse *fat)
{
	// Traverses until first position of next line
	// and eventually tries to refill array

	daStream *das = fat->das;
	while(LFCR[ (unsigned)*fat->das->r_iter ] == zval)
	{
		++das->r_iter;
		if(fatCheckFill(fat))
			return fat->state;
	}
	//printf("[fat_skipLine] post iter: '%s'\n",das->r_iter);
	return fat_skipNewLine(fat);
}

static R_INLINE int fat_skipComment(faTraverse *fat)
{
	if(*fat->das->r_iter == fat_char_cmt)
	{
		++fat->nCom;
		if(fat_skipLine(fat))
		{
			fat->state &= (~fat_loc);
			return fat_ok;
		}
	}
	return fat->state;
}

static R_INLINE int fat_skipN(faTraverse *fat)
{
	// Proceeds in array as long as 'N' nucleotide is present.
	// Tries to skip Newline (LF) characters

	// Eventually checks for new sequence and eventually
	// sets flag (No initializing on new sequence).

	// Eventually tries to refill array.

	if(fatEmpty(fat))
		return fat->state;

	daStream *das=fat->das;
	while(ACGTN[ (unsigned)*das->r_iter ] == nval)
	{
		++das->r_iter;
		++fat->nN;

		if(fatCheckFill(fat))
			return fat->state;

		//printf("[skip  N] Post fill '%s'\n",das->r_iter);

		// Function does not return
		// with r_iter on newline
		if(fat_checkNewLine(fat))
			if(fat_skipNewLine(fat))
				return fat->state;
	}

	// Unset N-flag
	fat->state &= (~fat_loc);
	return fat_ok;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Routines:
// 			fat_getSeqName
//			fat_findNextKmer
//			fat_initSeq
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

char * fat_getSeqName(faTraverse *fat)
{
	// From actual sequence delimiter:
	// Copies text until empty space (or newline)
	// into new char array.

	// Char array will be returned und must be free'd from outside.

	// Eventually tries to refill array
	daStream *das = fat->das;
	unsigned nchar = 0;

	//printf("[fat_getSeqName] '%s'\n",fat->das->r_iter);

	if(*das->r_iter == fat_char_newSeq)
	{
		char *seq_name = das->r_iter;
		while( !( (*seq_name==' ') || (*seq_name==fat_char_lf) ) )
		{
			++seq_name;
			++nchar;

			// End of array reached
			//printf("[fat_getSeqName]  seq_name: '%s'\n",seq_name);
			if(seq_name == das->r_end)
			{
				if(das_fill(das))
				{
					//printf("[fat_getSeqName] empty seq name.\n");
					return 0;
				}
				// Restart search
				seq_name = fat->das->r_iter;
				nchar = 0;
			}
		}

		if(nchar == 0)
		{
			//printf("[fat_getSeqName] empty seq name.\n");
			return 0;
		}

		// Skip first char='>'
		--nchar;
		seq_name=fat->das->r_iter;
		++seq_name;

		//printf("[fat_getSeqName] copy: '%s'\tnchar=%u\n",das->r_iter,nchar);

		// Must be free'd from outside
		char *ret = calloc(sizeof(char), nchar + 1);
		unsigned i;

		for(i=0; i<nchar; ++i)
		{
			ret[i] = *seq_name;
			++seq_name;
		}
		//printf("[seq name] success.\t'%s'\n",fat->das->r_iter);
		return ret;
	}
	return 0;
}


// ToDo: Change into loop entry point from top_level
int fat_skipSeqHeader(faTraverse *fat)
{
	if(fat->state & fat_newSeq)
	{
		++fat->nSeq;
		fat->state &= (~fat_newSeq);
		// Skip header line
		if(fat_skipLine(fat))
			return fat->state;
	}
	return fat_ok;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Routines:
// 			fat_returnFromTranfer
//			fat_procNuc
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE int fat_retProcNuc(faTraverse *fat)
{
	// Joint terminating routine for
	// procNuc function.
	daStream *das=fat->das;
	//printf("[ret proc] plen: %lu\t'%s'\n",(das->p_iter-das->pos),das->pos);

	if( (das->p_iter-das->pos) >= (fat->k) )
	{
		// A)	At least k nucleotides copied
		// 		-> Success:
		//		Add string terminator.
		if((das->p_iter-das->pos) > 0)
		{
			das->npPos = (int)(das->p_iter - das->pos);
			*das->p_iter = fat_char_eos;
		}

		//		Set nuc_ready state.
		fat->state |= fat_nucReady;
		fat->state &= (~fat_loc);

		//printf("[proc nuc] NUC READY: '%s'\n",das->pos);
		return fat->state;
	}
	else
	{
		// B)	Not enough nucleotides
		//		four counting
		//		-> Failure:
		//		Reset pos array.
		das->p_iter = das->pos;
		*das->p_iter = '\0';
		//		Set pos empty state
		fat->state &= (~fat_nucReady);
		fat->state |= fat_loc;
		//printf("[proc nuc] Nuc_empty\n");
		return fat_ok;
	}
}


int fat_procNuc(faTraverse *fat)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Copies continuous Nucleotide sequence from
	// input  fasta stream		(rfc array) to
	// output nucleotide stream	(pos array).
	//
	// Tries to recover when newline character is found or
	// end of rfc-array is reached.
	//
	// Copying always starts at begin of nuc array
	//
	// Returns fat_ok (=0) when *NO* output is created.
	//
	// Allows 'if(fat_procNuc) => process output'
	// testing.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	daStream *das = fat->das;
	fat->state &= (~fat_nucReady);

	//printf("[proc Nuc] start: '%s'\n",das->r_iter);
	if((ACGT[ (unsigned)*das->r_iter ] != zval))
	{
		// Goto begin of proc array.
		das->p_iter = das->pos;

		//unsigned i=0,j;
		while(!dasEmpty(das))
		{
			// Copy nucleotides
			while( !( (ACGT[ (unsigned)*das->r_iter ] == zval) || dasEmpty(das) || dasProcEmpty(das) ) )
			{
				*das->p_iter = *das->r_iter;
				++das->r_iter;
				++das->p_iter;
			}

			// Two recoverable conditions:
			if(fat_checkNewLine(fat))
			{
				//printf("[proc Nuc] New Line.\n");
				if(fat_skipNewLine(fat))
					return fat_retProcNuc(fat);
			}
			else if(fatEmpty(fat)) // Can't be checked with fatCheckFill
			{
				//printf("[proc Nuc] fatEmpty.\n");
				if(das_fill(das))
					return fat_retProcNuc(fat);
			}
			else // Anything else -> break
			{
				//printf("[proc Nuc] else ...\n");
				return fat_retProcNuc(fat);
			}
			//printf("[proc Nuc] End while pos: '%s'\n",das->pos);
		}

		// Pos array is full
		fat->state |= fat_nucReady;
		fat->state &= (~fat_loc);
		//printf("[proc Nuc] return fat_ok\n");
		return fat_retProcNuc(fat);
		//return fat->state;
	}
	else //if(fat_checkNewLine(fat))
	{
		//printf("[proc Nuc] AGCT->zval\n");
		if(fat_checkNewLine(fat))
			fat_skipNewLine(fat);

			/*
			if(fat_skipNewLine(fat))
			{
				Rprintf("[proc Nuc] Newline skipped: '%s'\n", das->r_iter);
			}
			else
				Rprintf("[proc Nuc] Found ASCII %u\n", (unsigned char)(*fat->das->r_iter));
		*/
	}

	fat->state |= fat_loc;
	return fat->state;
}



int fat_findKarray(faTraverse *fat)
{
	//fat->state&=(~fat_nucReady);
	//printf("[fat_findKarray] Start: '%s'\t%u\tstrlen: %lu\n",fat->das->r_iter,fatEmpty(fat),strlen(fat->das->r_iter));
	if(!fatEmpty(fat))
	{
		if(fat_checkN(fat))
		{
			//printf("[fat_findKarray] skipN.\n");
			fat_skipN(fat);
		}
		if(fat_checkNewSeq(fat))
		{
			//printf("[fat_findKarray] New seq found.\n");
			return fat->state;
		}
		if(fat_checkComment(fat))
		{
			fat_skipComment(fat);
			//printf("[fat_findKarray] post cmt: '%s'\n",fat->das->r_iter);
		}
		if(fat_procNuc(fat))
		{
			//printf("[fat_findKarray] procNuc return state: %i\n",fat->state);
			return fat->state;
		}
		/*
		if(fat_checkNewLine(fat))
		{
			printf("[fat_findKarray] CheckNewLine!\n");
		}
		*/
	}
	//printf("[find Kar] Something else: '%s'\n",fat->das->r_iter);
	return fat_ok;
}



#endif /* FA_TRAVERSE_H_ */
