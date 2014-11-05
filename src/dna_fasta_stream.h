/*
 * dna_fasta_stream.h
 *
 *  Created on: 02.10.2013
 *      Author: wolfgang
 */

#ifndef DNA_FASTA_STREAM_H_
#define DNA_FASTA_STREAM_H_

# include <zlib.h>
# include "stat_defs.h"

// Assumes that line size >= k

/*
 *
 * http://en.wikipedia.org/wiki/FASTA_format
 *
 * The description line is distinguished from the sequence data by a greater-than (">")
 * symbol in the first column.
 * The word following the ">" symbol is the identifier of the sequence
 * and the rest of the line is the description (both are optional).
 * There should be no space between the ">" and the first letter of the identifier.
 * It is recommended that all lines of text be shorter than 80 characters.
 */



///////////////////////////////////////////////////////////////////////////////////////////////////
//
// File stream (can handle text and gz files)
//
///////////////////////////////////////////////////////////////////////////////////////////////////


static const unsigned faf_file_open		= 0;
static const unsigned faf_file_eof		= 1;
static const unsigned faf_file_closed	= 2;
static const unsigned faf_txt    		= 4;
static const unsigned faf_gz     		= 8;
static const unsigned faf_error  		= 16;

typedef struct fasta_file_stream
{
	unsigned type;
	unsigned status;

	FILE *f;
	gzFile gz;

} fafStream;

static inline unsigned faf_isEof(fafStream *faf) { return faf->status & faf_file_eof; }
static inline unsigned faf_isOpen(fafStream *faf)
{
	if(faf->status == faf_file_open)
		return 1;
	return 0;
}

static fafStream* faf_stream_init(const char* filename, unsigned mode)
{
	// Construct fafStream object from opened file
	fafStream *faf = calloc(sizeof(fafStream), 1);
	if(mode == faf_gz)
	{
		faf->type = faf_gz;
		faf->gz = gzopen(filename,"rb");
		if(faf->gz)
			faf->status = faf_file_open;
		else
			faf->status = faf_file_closed;
	}
	else
	{
		faf->type = faf_txt;
		faf->f = fopen(filename, "r");
		if(faf->f)
			faf->status = faf_file_open;
		else
			faf->status = faf_file_closed;
	}
	return faf;
}

static void faf_destroy(fafStream *faf)
{
	if(!(faf->status & faf_file_closed))
	{
		if(faf->type == faf_gz)
			gzclose(faf->gz);
		else
			fclose(faf->f);
	}
	free(faf);
}

static size_t inline faf_read(fafStream *faf, char *dest, unsigned size)
{
	if(faf->status == faf_file_open)
	{
		unsigned nchar;
		if(faf->type == faf_gz)
		{
			nchar = gzread(faf->gz, dest, sizeof(char) * size);

			if(gzeof(faf->gz))
			{
				faf->status = faf_file_eof;
				gzclose(faf->gz);
				faf->gz = 0;
				faf->status |= faf_file_closed;
			}
		}
		else
		{
			nchar = fread(dest, sizeof(char), size, faf->f);
			if(feof(faf->f))
			{
				faf->status = faf_file_eof;
				fclose(faf->f);
				faf->f = 0;
				faf->status |= faf_file_closed;
			}
		}
		return nchar;
	}
	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// fasta Stream
//
///////////////////////////////////////////////////////////////////////////////////////////////////


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Define buffer size and delimiting characters
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static const unsigned long	fas_loc_Nchar	= 10;
static const unsigned long	fas_size		= 10+1;
static const char			fas_seq_delim	= '>';
static const char			fas_comment		= ';';
static const char			fas_loc_NewLine	= '\n';
static const char			fas_eoc			= '\0';


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Define status flags
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static const int fas_err				= -1;	// Unrecoverable error state
static const int fas_ok					= 0;	// Must be 0 !!!

static const int fas_loc_kReady			= 1;
static const int fas_nuc_ready          = 2;
// Recoverable reasons for transfer interruption
static const int fas_loc_newLine		= 4;
static const int fas_loc_newSeq			= 8;
static const int fas_loc_comment		= 16;
static const int fas_loc_N				= 32;
// File and stream running out of content
static const int fas_stream_eof			= 64;
static const int fas_stream_empty	    = 128;
static const int fas_nuc_empty		    = 256;

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// 	Struct definition
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

typedef struct fasta_stream
{

	fafStream *fasta;

	char * fas;			// Raw fasta content
	char * nuc;			// Processed pure Nucleotide sequence

	char * iter;		// iterator for fas
	char * nuc_iter;	// iterator for nuc

	unsigned ffc;		// future fas content = fas_loc_Nchar - (iter-fas)
	unsigned pfc;		// past   fas content = iter-fas (on demand)

	unsigned fnc;		// future nuc content
	unsigned pnc;		// past   nuc content

	int stream_state;	// Carries stream state flags

	unsigned nN;		// Number of N's so far seen
	unsigned nSeq;		// Number of seq's so far seen
	unsigned k;			// K-mer length

} faStream;

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Constructing and and File operations
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

void faStream_destroy(faStream *fa)
{
	faf_destroy(fa->fasta);
	fa->fasta = 0;

	free(fa->fas);
	fa->fas = 0;

	free(fa->nuc);
	fa->nuc = 0;

	free(fa);
}


int fas_fill(faStream *fa)
{
	if(!faf_isEof(fa->fasta))
	{
		// Unused suffix in array?
		if(fa->ffc > 0)
		{
			fa->pfc = fa->iter - fa->fas;

			// Enough space at begin of array?
			if(fa->pfc < fa->ffc)
			{
				fa->stream_state = fas_err;
				return fa->stream_state;
			}

			// Shift unused suffix to begin
			memcpy(fa->fas, fa->iter, fa->ffc);

			// Fill rest of array from file
			fa->ffc += faf_read(fa->fasta, fa->fas + fa->ffc, fa->pfc);
		}
		else
		{
			// Refill whole array
			fa->ffc = faf_read(fa->fasta, fa->fas, fas_loc_Nchar);
		}
		fa->iter = fa->fas;
	}

	// Set flags
	if(faf_isEof(fa->fasta))
		fa->stream_state |= fas_stream_eof;

	if(fa->ffc == 0)
	{
		fa->stream_state |= fas_stream_empty;
		return fa->stream_state;
	}

	// Return success
	fa->stream_state &= (~fas_stream_empty);
	return fas_ok;
}

faStream * faStream_init(const char* filename, unsigned k, unsigned mode)
{
	if(k > fas_size)
	{
		printf("[faStream_init] k > fas_size!\n");
		return 0;
	}

	faStream *fa = calloc(sizeof(faStream), 1);

	if(mode == faf_gz)
		fa->fasta = faf_stream_init(filename, faf_gz);
	else
		fa->fasta = faf_stream_init(filename, faf_txt);

	if(!faf_isOpen(fa->fasta))
	{
		printf("[faStream_init] Cannot open file '%s'!\n", filename);
		faStream_destroy(fa);
		return 0;
	}

	fa->k=k;
	fa->fas = calloc(fas_size, sizeof(char));
	fa->nuc = calloc(fas_size, sizeof(char));

	// fas_fill will read whole array
	fa->stream_state = fas_stream_empty;
	fa->ffc = 0;

	// Provide initial filling
	if(!fas_fill(fa))
	{
		printf("[faStream_init] Initial array filling failed!\n");
		faStream_destroy(fa);
		return 0;
	}
	return fa;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Check routines
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static inline int fas_fas_end(faStream *fa)		{ return fa->ffc == 0; }
static inline int fas_nuc_end(faStream *fa)		{ return fa->fnc == 0; }



// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Skipping routines for Line feed, entire line and 'N' Nucleotides
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


static inline int fas_checkNewLine(faStream *fa)
{
	if(*fa->iter == fas_loc_NewLine)
	{
		fa->stream_state |= fas_loc_newLine;
		return fa->stream_state;
	}
	return fas_ok;
}

static inline int fas_skipNewLine(faStream *fa)
{
	// Skips LF ('\n') character and
	// eventually tries to fill array

	if(*fa->iter == fas_loc_NewLine)
	{
		if(fa->ffc == 0)
		{
			if(!fas_fill(fa))
				return fa->stream_state;
		}
		// Eat newline
		++fa->iter;
		--fa->ffc;
		// Unset flag
		fa->stream_state &= (~fas_loc_newLine);
		return fas_ok;
	}
	return fa->stream_state;
}


static inline int fas_skipLine(faStream *fa)
{
	// Traverses until first position of next line
	// and eventually tries to refill array
	while(*fa->iter != fas_loc_NewLine)
	{
		if(fa->ffc > 0)
		{
			++fa->iter;
			--fa->ffc;
		}
		else
		{
			if(!fas_fill(fa))
				return fa->stream_state;
		}
	}
	return fas_skipNewLine(fa);
}

static inline int fas_checkNewSeq(faStream *fa)
{
	// Checks for Seq-delimiter at current position
	// and eventually sets flag
	// (does not proceed)
	if(*fa->iter == fas_seq_delim)
	{
		fa->stream_state |= fas_loc_newSeq;
		return fa->stream_state;
	}
	return fas_ok;
}

static inline int fas_skipSeqHeader(faStream *fa)
{
	if(*fa->iter==fas_seq_delim)
	{
		if(!fas_skipLine(fa))
		{
			fa->stream_state&=(~fas_loc_newSeq);
			return fas_ok;
		}
	}
	return fa->stream_state;
}

static inline int fas_checkComment(faStream *fa)
{
	if(*fa->iter == fas_comment)
	{
		fa->stream_state |= fas_loc_comment;
		return fa->stream_state;
	}
	return fas_ok;
}

static inline int fas_skipComment(faStream *fa)
{
	if(*fa->iter == fas_comment)
	{
		if(!fas_skipLine(fa))
		{
			fa->stream_state &= (~fas_loc_comment);
			return fas_ok;
		}
	}
	return fa->stream_state;
}

static inline int fas_checkN(faStream *fa)
{
	if(ACGTN[ (unsigned)*fa->iter ] == nval)
	{
		fa->stream_state |= fas_loc_N;
		return fa->stream_state;
	}
	return fas_ok;
}

static inline int fas_skipN(faStream *fa)
{
	// Proceeds in array as long as 'N' nucleotide is present.
	// Tries to skip Newline (LF) characters

	// Eventually checks for new sequence and eventually
	// sets flag (No initializing on new sequence).

	// Eventually tries to refill array.

	while(ACGTN[(unsigned)*fa->iter] == nval)
	{
		++fa->iter;
		--fa->ffc;
		++fa->nN;

		// Array exhausted
		if(fa->ffc == 0)
		{
			if(!fas_fill(fa))
				return fa->stream_state;
		}
		// Newline found
		if(!fas_checkNewLine(fa))
		{
			if(!fas_skipNewLine(fa))
				return fa->stream_state;

			while(!fas_checkComment(fa))
			{
				if(!fas_skipComment(fa))
					return fa->stream_state;
			}
			if(!fas_checkNewSeq(fa))
				return fa->stream_state;
		}
	}
	// Unset N-flag
	fa->stream_state &= (~fas_loc_N);
	return fas_ok;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Routines:
// 			fas_getSeqName
//			fas_findNextKmer
//			fas_initSeq
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

char * fas_getSeqName(faStream *fa)
{
	// From actual sequence delimiter:
	// Copies text until empty space (or newline)
	// into new char array.

	// Char array will be returned und must be free'd from outside.

	// Eventually tries to refill array

	if(*fa->iter == fas_seq_delim)
	{
		char *seq_name = fa->iter;
		unsigned ffc = fa->ffc;

		while(!( (*seq_name==' ') || (*seq_name==fas_loc_NewLine) || (ffc==0) ) )
		{
			++seq_name;
			--ffc;

			// End of array reached
			if(ffc==0)
			{
				// Try to refill
				if(!(fa->stream_state & fas_stream_eof))
					fas_fill(fa);
				else
					return 0;

				// Restart search
				seq_name = fa->iter;
				ffc = fa->ffc;
			}
		}

		// Skip first char='>'
		unsigned i, nchar;

		nchar = seq_name-fa->iter - 1;

		// Must be free'd from outside
		char *ret = calloc(sizeof(char), nchar + 1);

		seq_name = fa->iter;
		++seq_name;
		for(i=0; i<nchar; ++i)
		{
			ret[i] = *seq_name;
			++seq_name;
		}
		return ret;
	}
	return 0;
}

int fas_initSeq(faStream *fa)
{
	// + + + + + + + + + + + + + + + + + + + + + //
	// Performs initializing steps when new
	// sequence delimiter ('>') is found:
	//
	// A) Clears newSeq flag
	// B) Skip actual line
	// C) Skip following comment lines
	// D) Calls findNextKmer
	// (there is no testing for comments
	// at other positions).
	// + + + + + + + + + + + + + + + + + + + + + //

	// Clear newSeq flag
	fa->stream_state &= (~fas_loc_newSeq);

	if(*fa->iter == fas_seq_delim)
	{
		// Skip header line
		if(fas_skipLine(fa) != fas_ok)
			return fa->stream_state;

		// Skip comment lines
		while(*fa->iter == fas_comment)
		{
			if(fas_skipLine(fa) != fas_ok)
				return fa->stream_state;
		}
	}
	return fas_ok;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Routines:
// 			fas_returnFromTranfer
//			fas_TransferNucArray
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static inline int fas_returnFromTranfer(faStream *fa)
{
	// Joint terminating routine for
	// TransferNucArray function.

	if((fa->nuc_iter - fa->nuc) >= fa->k)
	{
		// A)	At least k nucleotides copied
		// 		-> Success:
		//		Add string terminator.
		if(fa->fnc > 0)
			*fa->nuc_iter = fas_eoc;

		//		Set nuc_ready state.
		fa->stream_state |= fas_nuc_ready;
		fa->stream_state &= (~fas_nuc_empty);
		return fas_ok;
	}
	else
	{
		// B)	Not enough nucleotides
		//		four counting
		//		-> Failure:
		//		Reset nuc array.
		fa->nuc_iter = fa->nuc;
		*fa->nuc = '\0';

		//		Set nuc empty state
		fa->stream_state &= (~fas_nuc_ready);
		fa->stream_state |= fas_nuc_empty;
		return fa->stream_state;
	}
}

int fas_TransferNucArray(faStream *fa)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Copies continuous Nucleotide sequence from
	// input  fasta stream		(fas array) to
	// output nucleotide stream	(nuc array).
	//
	// Tries to recover when newline character is found or
	// end of fas-array is reached.
	//
	// Transaction skips newline characters and terminates
	// when either 'N' or '>' (sequence delimiter)
	// is found.
	//
	// Copying always starts at begin of nuc array
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	// Goto begin of nuc array.
	fa->nuc_iter = fa->nuc;
	fa->fnc = fas_loc_Nchar;


	while(!fas_nuc_end(fa))
	{
		// Copy nucleotides
		while( (ACGT[ (unsigned)*fa->iter ] != zval) && (fa->ffc > 0) && (fa->fnc > 0) )
		{
			*fa->nuc_iter = *fa->iter;
			++fa->iter;
			++fa->nuc_iter;
			--fa->ffc;
			--fa->fnc;
		}

		if(fas_checkNewLine(fa))
		{
			if(!fas_skipNewLine(fa))
				return fas_returnFromTranfer(fa);

			if(fas_checkNewSeq(fa))
				return fas_returnFromTranfer(fa);

			while(fas_checkComment(fa))
			{
				if(fas_skipComment(fa))
					return fas_returnFromTranfer(fa);
			}
			// Here: NewLine readily skipped
		}
		else if(fas_fas_end(fa))
		{
			if(!fas_fill(fa))
				return fas_returnFromTranfer(fa);
		}
		else if(fas_checkN(fa))
			return fas_returnFromTranfer(fa);
	}

	// Nuc array is full
	fa->stream_state |= fas_nuc_ready;
	return fas_ok;
}


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Routines: Accessors for stream flags.
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


int fa_empty(faStream *fa)    { return (fa->stream_state & fas_stream_empty); }
int fa_NucReady(faStream *fa) { return (fa->stream_state & fas_nuc_ready);    }
int fa_K_Ready(faStream *fa)  { return (fa->stream_state & fas_loc_kReady);   }
int fa_NewSeq(faStream *fa)   { return (fa->stream_state & fas_loc_newSeq);   }
int fa_N_Nuc(faStream *fa)    { return (ACGTN[(unsigned)*fa->iter] == nval);  }

// Is unset when Nuc-Array is processed.
void fa_unsetNucReady(faStream *fa)
{
	fa->stream_state &= (~fas_nuc_ready);
	fa->stream_state |= fas_nuc_empty;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


#endif /* DNA_FASTA_STREAM_H_ */
