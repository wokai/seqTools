/*
 * fastq_parser.h
 *
 *  Created on: 20.10.2013
 *      Author: wolfgang
 */

#ifndef FASTQ_PARSER_H_
#define FASTQ_PARSER_H_

# include <zlib.h>

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Define buffer size and delimiting characters
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
static const unsigned	fqp_gz_size			= 0xFFFFFF;	// 16.777.215 gzFile
													// internal buffer size

static const char	fqp_char_sqDelim 	= '@';		// FASTQ sequence delimiter
static const char	fqp_char_qlDelim 	= '+';		// FASTQ quality  delimiter
static const char	fqp_char_lf			= '\n';		// Line feed
static const char	fqp_char_eos		= '\0';		// End of string


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Define file status flags
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
static const int	fqp_file_closed = 0;
static const int	fqp_file_open	= 1;


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// Define FASTQ format status flags
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
static const int 	fqp_fastq_ok	= 0;
static const int	fqp_fastq_error	= -1;


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 	Struct definition
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

typedef struct fastq_file_parser
{

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	// File access		declarations
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	gzFile gz;
	int fstate;	// File state

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	// Character array	declarations
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	int cap;		// Capacity: Maximal number of characters in rbuf and pbuf
						// (array size=cap+1)

	// Pointers to character arrays (= begin pointers)
	char * rbuf;		// Buffer for raw file content
	char * pbuf;		// Buffer for processed output sequence

	// Iterators which reside inside arrays
	char * r_iter;		// iterator for rbuf
	char * p_iter;		// iterator for pbuf

	//  Point to past-the-end character of arrays ('\0')
	char * r_end;		// rbuf array
	char * p_end;		// behind last retrievable character


	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	// Array content	declarations
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  //
	// Number of character values in pbuf array (=strlen)
	int seqlen;

	int fastq_format_state;

	// Internal counters
	int nSeq;
	unsigned minSeqLen;
	unsigned maxSeqLen;
	int nN;


} fqParser;

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Check routines
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE int fqpEmpty(fqParser *fqp)		    { return (fqp->r_iter == fqp->r_end); }
static R_INLINE int fqpProcFull(fqParser *fqp)	    { return (fqp->p_iter == fqp->p_end); }
static R_INLINE int fqpIsOpen(fqParser *fqp)		{ return (fqp->fstate == fqp_file_open); }
static R_INLINE int fqpFastqError(fqParser *fqp)
								{ return (fqp->fastq_format_state == fqp_fastq_error); }

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Constructing and and File operations
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE int fqp_open(fqParser *fqp,const char *filename)
{
	fqp->gz = gzopen(filename, "rb");

	if(!fqp->gz)
	{
		fqp->fstate = fqp_file_closed;
		return 0;
	}


	/*
	// Increase buffer size (from 8192 bytes)
	if(gzbuffer(fqp->gz,fqp_gz_size)<0)
	{
		free(fqp);
		fqp->gz=0;
		fqp->fstate=fqp_file_closed;
		return 0;
	}
	*/

	fqp->fstate = fqp_file_open;
	return 1;
}

static R_INLINE void fqp_close(fqParser *fqp)
{
	if(fqp->fstate == fqp_file_open)
	{
		if(fqp->gz)
		{
			gzclose(fqp->gz);
			fqp->gz = 0;
		}
		fqp->fstate = fqp_file_closed;
	}
}

static R_INLINE void fqp_print_error(fqParser *fqp)
{
	if(fqp->fstate == fqp_file_open)
	{
		if(fqp->gz)
		{
			int err;
			const char * const msg=gzerror(fqp->gz, &err);
			if(err)
				Rprintf("[FASTQ file ERROR] GZ: %s\n",msg);
		}
		else
			Rprintf("[FASTQ file ERROR] File pointer=0.\n");
	}
	else
		Rprintf("[FASTQ file ERROR] File closed.\n");

}

void fqp_destroy(fqParser *fqp)
{
	if(fqp)
	{
		if(fqp->fstate == fqp_file_open)
			fqp_close(fqp);

		free(fqp->rbuf);
		fqp->rbuf = 0;

		free(fqp->pbuf);
		fqp->pbuf = 0;

		free(fqp);
	}
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Basic buffer operations
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE void fqp_initRbuf_ptrs(fqParser *fqp)
{
	// last array position, contains '\0'
	fqp->r_end = fqp->rbuf+fqp->cap;

	// Entails complete filling
	fqp->r_iter = fqp->r_end;
}

static R_INLINE void fqp_initPbuf_ptrs(fqParser *fqp)
{
	fqp->p_end = fqp->pbuf + fqp->cap;
	fqp->p_iter = fqp->pbuf;
}

static R_INLINE void fqp_clearPbuf(fqParser *fqp)
{
	fqp->p_iter = fqp->pbuf;
	*fqp->p_iter = fqp_char_eos;
}

static R_INLINE void fqp_clearRbuf(fqParser *fqp)
{ fqp->r_iter = fqp->r_end; }

int fqp_getR_content_size(fqParser *fqp)
{ return (int) (fqp->r_end - fqp->r_iter); }


int fqp_fill_rbuf(fqParser *fqp)
{
	int nRead;
	unsigned nDemand;

	if(fqp->fstate == fqp_file_open)
	{
		if(fqp->r_iter < fqp->r_end)
		{
			if(fqp->r_iter > fqp->rbuf)
			{
				char* rhs = fqp->r_iter;
				fqp->r_iter = fqp->rbuf;

				while(rhs != fqp->r_end)
				{
					*fqp->r_iter = *rhs;
					++fqp->r_iter;
					++rhs;
				}
			}
		}
		else
		{
			fqp->r_iter = fqp->rbuf;
		}

		nDemand = (unsigned)(fqp->r_end - fqp->r_iter);
		nRead = gzread(fqp->gz, fqp->r_iter, sizeof(char) * nDemand);
		if(nRead < ((int) nDemand) )
		{
			if(!gzeof(fqp->gz))
				fqp_print_error(fqp);

			fqp_close(fqp);

			// Always left shift
			fqp->r_end = fqp->r_iter + nRead;

			// Place '\0'
			*fqp->r_end = fqp_char_eos;
		}
		// After Filling iter must point to begin of array
		fqp->r_iter = fqp->rbuf;

		return nRead;
	}
	return 0;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// Initialize structure
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

fqParser * fqp_init(const char* filename, int buf_capacity)
{

	// + + + + + + + + + + + + + + + + + + + + + //
	// Initialize raw structure memory
	// + + + + + + + + + + + + + + + + + + + + + //
	fqParser *fqp = calloc(sizeof(fqParser), 1);
	if(!fqp)
		return 0;

	// + + + + + + + + + + + + + + + + + + + + + //
	// File open
	// + + + + + + + + + + + + + + + + + + + + + //
	if(!fqp_open(fqp, filename))
			return fqp;

	// + + + + + + + + + + + + + + + + + + + + + //
	// Initialize buffer memory
	// + + + + + + + + + + + + + + + + + + + + + //
	fqp->cap = buf_capacity; // fqp_buf_size
	size_t buf_size = ((size_t) fqp->cap) + 1;

	fqp->rbuf = calloc(buf_size, sizeof(char));
	if(!fqp->rbuf)
	{
		fqp_destroy(fqp);
		return 0;
	}
	fqp->pbuf = calloc(buf_size, sizeof(char));
	if(!fqp->pbuf)
	{
		fqp_destroy(fqp);
		return 0;
	}

	// + + + + + + + + + + + + + + + + + + + + + //
	// Initialize buffer pointer
	// + + + + + + + + + + + + + + + + + + + + + //
	fqp_initRbuf_ptrs(fqp);
	fqp_initPbuf_ptrs(fqp);

	// + + + + + + + + + + + + + + + + + + + + + //
	// Initialize counter values
	--fqp->minSeqLen;

	//fqp_fill_rbuf(fqp);

	// + + + + + + + + + + + + + + + + + + + + + //
	// Ready:
	// - Returns memory initialized structure
	// - File is possibly not open
	//   (file not found)
	// - When file is open, first buffer filling
	//   is done.
	// + + + + + + + + + + + + + + + + + + + + + //
	return fqp;
}


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
//	Fill operations
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

int fqpCheckFill(fqParser *fqp)
{
	if(fqp->r_end == fqp->r_iter)
		return fqp_fill_rbuf(fqp);
	return 0;
}

int fqpCheckCapFill(fqParser *fqp, unsigned capacity)
{
	if( (fqp->r_end - fqp->r_iter) < capacity)
		return fqp_fill_rbuf(fqp);
	return 0;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
//	Skip routines
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE int fqp_checkSkipNewLine(fqParser *fqp)
{
	if(*fqp->r_iter == fqp_char_lf)
	{
		fqpCheckCapFill(fqp, 2);

		if(fqp->r_iter < fqp->r_end)
			++fqp->r_iter;
		return 1;
		// if(skip) => newLine skipped
	}
	return 0;
}

static R_INLINE int fqp_skipLine(fqParser *fqp)
{
	// Traverses until first position of next line
	// and eventually tries to refill array
	while( (*fqp->r_iter != fqp_char_lf) && (fqp->r_iter < fqp->r_end) )
	{
		++fqp->r_iter;
		fqpCheckFill(fqp);
	}
	if(fqp->r_iter<fqp->r_end)
	{
		++fqp->r_iter;
		fqpCheckFill(fqp);
		return 1;
	}
	// if (skip) => newLine skipped
	return 0;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
//	Copy routines
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE void fqp_copyLine(fqParser *fqp)
{
	// Does not reset p_iter !
	while( !( (*fqp->r_iter == fqp_char_lf) || (fqp->r_end == fqp->r_iter) || (fqp->p_iter == fqp->p_end) ) )
	{
		*fqp->p_iter = *fqp->r_iter;
		++fqp->r_iter;
		++fqp->p_iter;
		fqpCheckFill(fqp);
	}
	if( (*fqp->r_iter == fqp_char_lf) && (fqp->r_iter != fqp->r_end))
	{
		++fqp->r_iter;
		fqpCheckFill(fqp);
	}
	fqp->seqlen = (int)(fqp->p_iter-fqp->pbuf);
}

static R_INLINE int fqp_copyNchar(fqParser *fqp, int nchar)
{
	// RESET'S p_iter !
	fqp->p_iter = fqp->pbuf;

	if(nchar > fqp->cap)
		return 0;

	int nDemand = nchar-(int)(fqp->r_end - fqp->r_iter);

	if(nDemand > 0)
	{
		if(fqp_fill_rbuf(fqp) < nDemand)
		{
			fqp_clearRbuf(fqp);
			fqp_clearPbuf(fqp);
			return 0;
		}
	}

	int nCopied = 0;
	while( !( (nCopied == nchar) || (fqp->r_iter == fqp->r_end) ) )
	{
		if(*fqp->r_iter == fqp_char_lf)
		{
			++fqp->r_iter; // skip newline
			fqpCheckFill(fqp);
		}
		else
		{
			*fqp->p_iter = *fqp->r_iter;
			++fqp->p_iter;
			++fqp->r_iter;
			++nCopied;
			fqpCheckFill(fqp);
		}
	}
	return nCopied;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
//	Process Nucleotide sequence and Phred Quality sequence
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


// Find End of line
// Used for locating the end of the header line
// in trim_fastq
size_t find_eol(const char * c)
{
	const char *iter = c;

	while(*iter != '\0' && *iter != '\n')
		++iter;

	return ( (size_t) (iter - c) );
}


int fqp_procNuc(fqParser *fqp)
{

	if(*fqp->r_iter != fqp_char_sqDelim)
	{
		//Rprintf("[proc nuc] fq format error.\n");
		fqp->fastq_format_state = fqp_fastq_error;
		return -1;
	}
	if(!fqp_skipLine(fqp))	// Header
		return 0;

	fqp_clearPbuf(fqp);
	while(!( (fqp->r_iter == fqp->r_end) || (*fqp->r_iter == fqp_char_qlDelim) || (fqp->p_iter == fqp->p_end) ) )
		fqp_copyLine(fqp);


	// + + + + + + + + + + + + + + + + + + + + + + //
	// Finalize copy procedure by adding
	// terminating '\0'
	// + + + + + + + + + + + + + + + + + + + + + + //
	if(fqp->p_iter != fqp->p_end)
	{
		if(fqp->p_iter != fqp->pbuf)
			++fqp->p_iter;

		*fqp->p_iter = fqp_char_eos;
	}
	// + + + + + + + + + + + + + + + + + + + + + + //

	if( (fqp->p_iter == fqp->p_end) && (*fqp->r_iter != fqp_char_qlDelim) )
	{
		Rprintf("\n[fastqq] Buffer overflow.\n");
		fqp->fastq_format_state = fqp_fastq_error;
		return -1;
	}

	if(fqp->seqlen > 0)
	{
		++fqp->nSeq;
	}

	return fqp->seqlen;
}

int fqp_procPhred(fqParser *fqp)
{
	if(*fqp->r_iter != fqp_char_qlDelim)
	{
		//Rprintf("[fastqq] Phred format error: '%s'\n",fqp->r_iter);
		fqp->fastq_format_state = fqp_fastq_error;
		return -1;
	}
	if(!fqp_skipLine(fqp))	// Header
		return 0;

	int nCopied = fqp_copyNchar(fqp, fqp->seqlen);
	// Assuming that the same number of chars are
	// copied as from fqp_procNuc the terminating
	// '\0' from before should be in right place.

	fqp_checkSkipNewLine(fqp);
	return nCopied;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


#endif /* FASTQ_PARSER_H_ */
