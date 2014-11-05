/*
 * dna_astream.h
 *
 *  Created on: 07.10.2013
 *      Author: wolfgang
 */

#ifndef DNA_ASTREAM_H_
#define DNA_ASTREAM_H_

#include "dna_fstream.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// dna file stream
//		Keeps two (equal sized) character buffers (raw and processed) and two iterators
//
//		Raw buffer can be (re-) filled from dna_fstream.
//		Client struct can define functions for searching and coypying content into
//		processed buffer.
//
///////////////////////////////////////////////////////////////////////////////////////////////////


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Define buffer size and delimiting characters
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static const char			das_char_eos			='\0';

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Define status flags
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static const int das_err		=-1;	// Unrecoverable error state
static const int das_ok			= 0;	// Must be 0 !!!
static const int das_empty		= 1;	// Future raw size =0

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// 	Struct definition
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

typedef struct dna_astream
{

	dfStream *dnaf;

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	// Character array declarations

	unsigned nchar;		// Maximal number of characters in rfc and pos
						// (array size=nchar+1)

	// Pointers to character arrays (= begin pointers)
	char * rfc;			// Raw file content				: rfc-array
	char * pos;			// Processed output sequence	: pos-array

	// Iterators which reside inside arrays
	char * r_iter;		// iterator for rfc
	char * p_iter;		// iterator for pos

	//  Point to past-the-end character of arrays ('\0')
	char * r_end;		// rfc array
	char * p_end;		// behind last retrievable character

	int npPos;		// Number of character values in pos array (=strlen)
	//
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //

	int state;	// Carries stream state flags

	unsigned nFill;
	unsigned nFillWhole;
	unsigned nFillWholeIncomp;
	unsigned nFillPart;
	unsigned nFillPartIncomp;

} daStream;

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Check routines
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static inline int dasEmpty(daStream *das) 		{ return das->r_iter == das->r_end; }
static inline int dasProcEmpty(daStream *das)	{ return das->p_iter == das->p_end; }
static inline int dasIsError(daStream *das)		{ return das->state & das_err; }
static inline int dasIsOpen(daStream *das)		{ return dfs_isOpen(das->dnaf); }
static inline int dasIsEof(daStream *das)		{ return dfs_stream_eof(das->dnaf); }

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Constructing and and File operations
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

void das_destroy(daStream *das)
{
	if(das)
	{
		dfs_destroy(das->dnaf);
		das->dnaf=0;
		free(das->rfc);
		das->rfc=0;
		free(das->pos);
		das->pos=0;
		free(das);
	}
}


int das_fill(daStream *das)
{
	//printf("[das_fill] before: '%s'\n",das->rfc);
	size_t count = 0;

	int rhs= (int)(das->r_end  - das->r_iter);		// Number of unprocessed characters in array
	int lhs= (int)(das->r_iter - das->rfc);			// Number of   processed characters in array

	if(dfs_isOpen(das->dnaf))
	{
		if(rhs>0)
		{
			//printf("[das_fill] rhs= %u\n",rhs);
			if(lhs<rhs)
			{
				// Enough space at begin of array?
				//printf("[das_fill] Partial fill ERROR!\n");
				das->state = das_err;
				return das->state;
			}
			// Shift unused suffix to begin
			memcpy(das->rfc, das->r_iter, (size_t) rhs);
			count = dfs_read(das->dnaf, das->rfc + rhs, (unsigned) lhs);

			/* lhs < 0 should not be possible */
			if(count < ((size_t) lhs))
				++das->nFillPartIncomp;
			++das->nFillPart;
		}
		else
		{	// Refill whole array
			//printf("[das_fill] rhs== %u\n",rhs);
			count = dfs_read(das->dnaf, das->rfc, das->nchar);
			if(count < das->nchar)
				++das->nFillWholeIncomp;
			++das->nFillWhole;
		}
		das->r_end = das->rfc + count; // past-the-end
		*das->r_end = das_char_eos;    // '\0'
		das->r_iter = das->rfc;        // Re-init iter
		++das->nFill;
	}

	if(count==0)
	{
		//printf("[das_fill] count==0.\n");
		das->state |= das_empty;
		return das->state;
	}
	// Return success
	das->state &= (~das_empty);
	return das_ok;
}

daStream * das_init(const char* filename, unsigned das_size)
{
	daStream *das = calloc(sizeof(daStream), 1);
	if(!das)
	{
		//printf("[das_init] das calloc returned 0!\n");
		return 0;
	}

	das->dnaf = dfs_stream_init(filename);

	if(!das->dnaf)
	{
		//printf("[das_init] dfs_stream_init returned 0!\n");
		das->state = das_err;
		return das;
	}

	das->nchar = das_size;
	das->rfc = calloc(das_size + 1, sizeof(char));
	if(!das->rfc)
	{
		//printf("[das_init] rfc calloc returned 0!\n");
		das->state = das_err;
		return das;
	}
	das->pos = calloc(das_size + 1, sizeof(char));
	if(!das->pos)
	{
		//printf("[das_init] pos calloc returned 0!\n");
		das->state = das_err;
		return das;
	}

	das->r_end = das->rfc + das_size;
	// Indicates empty buffer
	// -> first das_fill will
	// read complete buffer
	das->r_iter = das->r_end;

	// Returns memory initialized structure
	// but dfs file is possibly closed
	// (e.g. file not found).
	return das;
}


#endif /* DNA_ASTREAM_H_ */
