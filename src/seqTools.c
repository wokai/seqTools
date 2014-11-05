

#ifndef SEQTOOLS_C_
#define SEQTOOLS_C_
#include "seqTools.h"

/*
 * seqTools.c
 *
 *  Created on: 14.10.2013
 *      Author: kaisers
 */


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		create_dna_k_mers	:	Returns ordered vector with all DNA motifs for given k.
 * 								C-index of array is mindex.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


SEXP create_dna_k_mers(int len)
{
	int i, j, k, nstr;
	char * c;
	SEXP res;

	if(len > max_k)
		error("[create_dna_k_mers] Maximum value for pK is %u!", max_k);
	if(len <= 0)
		error("[create_dna_k_mers] pK must be > 0!");

	nstr= 1 << (len << 1);	/* =4^len */
	c = (char* ) R_alloc((size_t)(len + 1), sizeof(char));
	c[len] = '\0';

	res=PROTECT(allocVector(STRSXP, (R_xlen_t)nstr));
	for(i=0; i<nstr; ++i)
	{
		for(j=0, k=(2 * (len-1)); j<len; ++j)
		{
			c[j] = (char) rev_ACGT[(i >> k) & 3];
			k -= 2;
		}
		SET_STRING_ELT(res, (R_xlen_t)i, mkChar(c));
	}
	UNPROTECT(1);
	return res;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		do_countCheck_Kmers
 * 		k		:	Kmer length
 * 		inseq	:	Pointer into DNA sequence which to be counted
 * 		array	:	Array where motif counts are stored
 * 		npos	:	Number of positions = strlen(inseq) - k
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static R_INLINE bool do_countCheck_Kmers(char const * const inseq, int * const array, int * nn, int k, int npos)
{
	int 			m, j, i;
	unsigned		val_ok;
	unsigned long 	char_value, array_idx;
	char const 		* iter;

	/*
	 * Expected array length: 4^k
	 * iter traverses seq from * inseq for npos positions
	 */

	iter = inseq;
	m = k - 1;
	array_idx = 0;
	val_ok = 0;
	for(i=0; i<npos; ++i)
	{
		/*
		 * Calculate array_index for each Kmer at each seq-position
		 * and increase value at index-position by 1
		 */
		for(j=0; j<k; ++j)
		{
			/*
			 * Rprintf("[do_countCheck_Kmers] i=%u\tj=%i\tseq: %s\n", i, j, iter);
			 * Throw error when traversing '\0'
			 */
			if(iter[j] == '\0')
				error("[count_Kmer] Found string terminating NULL!");

			char_value = (unsigned long) ACGT[ (unsigned)iter[j] ];
			if(char_value == zval)
			{
				/*
				 *  Jump over 'N' character
				 */
				if(ACGTN[(unsigned)iter[j]] == nval)
				{
					++iter;
					array_idx = 0;
					++(*nn);

					/*
					 * current k-mer will not be counted
					 */
					val_ok = 0;
					break;
				}
				else
				{
					Rprintf("\n[do_countCheck_Kmers] Error : j: %u\t iter: '%s'\n", j, iter);
					return false;
				}
			}
			else
			{
				/*
				 * Index reverse order (because dna Kmers do it this way)
				 * (2* ((k-1)-j))= (m-j)<<1
				 */
				array_idx |= (char_value << ((m-j) << 1));
				/*
				 * current k-mer will be counted
				 */
				val_ok = 1;
			}
		}
		if(val_ok)
			++array[array_idx];
		++iter;
		array_idx = 0;
	}
	return true;
}

static R_INLINE int do_count_Kmers(char const * const inseq, int * const array, int *nn,int k, int npos)
{
	int i, j, m, np;
	unsigned char_value, array_idx, len, val_ok;
	char const * iter;

	/*
	 * Does not throw errors when Non-Nucs are found
	 * Advances for k positions when N is found
	 * k = length of Kmer
	 * if(k>max_k)
	 * 	error("[do_count_Kmers] k must be <= %u!", max_k);
	 */

	/*
	 * Expected array length: 4^k
	 * iter traverses seq from * inseq for npos positions
	 */

	len = (unsigned) strlen(inseq);

	np = npos;
	m = k - 1;
	if(len < (unsigned long)(np + m))
		return -1;

	iter = inseq;
	array_idx = 0;
	val_ok = 0;

	for(i=0; i<np; ++i)
	{

		/*
		 * Calculate array_index for each Kmer at each seq-position
		 * and increase value at index-position by 1
		 *
		 */

		/*
		 * Rprintf("[do_count_Kmers] Loop entry: i=%i\tseq: %s\n", i, iter);
		 */

		for(j=0; j<k; ++j)
		{
			/*
			 * Throw error when traversing '\0'
			 * if(iter[j]=='\0')
			 * return -1; //error("[do_count_Kmers] Found string terminating NULL!");
			 */

			char_value = ACGT[ (unsigned)iter[j] ];
			val_ok = 1;

			if(char_value == zval)
			{
				if(ACGTN[ (unsigned)iter[j] ] == nval)
				{
					iter += k;
					++(* nn);

					/* Past end of string?	*/
					if( ((unsigned)(iter-inseq)) >= len )
						return 0;

					/* k-1: compensate for already happened increase of i	*/
					np -= (k - 1);
					array_idx = 0;
					/* actual k-mer will not be counted */
					val_ok = 0;
				}
				break;
			}
			/* Index reverse order (because dna Kmers do it this way) */
			/* (2* ((k-1)-j))= (m-j)<<1 */
			array_idx |= (char_value << ((m-j) << 1));

			/*
			 * Rprintf("[do_count_Kmers] i=%i\tnp=%i\tj=%i\tidx=%u\tok=%u\tseq: %s\n", i, np, j, array_idx, val_ok, iter);
			 */

		}
		if(val_ok)
		{
			++array[array_idx];
			++iter;
		}
		array_idx = 0;
	}
	return 0;
}


static R_INLINE void count_first_Kmer(char const * const inseq, int * const array, int k)
{
	/*
	 * Counts one k-mer at beginning of inseq
	 * k = length of Kmer
	 */
	int j, m;
	unsigned long char_value, array_idx;

	m = k - 1;
	array_idx = 0;
	for(j=0; j<k; ++j)
	{
		/* Throw error when traversing '\0' */
		if(inseq[j] == '\0')
			error("[count_Kmer] Found string terminating NULL!");

		char_value = ACGT[ (unsigned)inseq[j] ];

		if(char_value == zval)
			return;
		/*
		 * Index reverse order (because dna Kmers do it this way)
		 * (2* ((k-1)-j))= (m-j)<<1
		 */
		array_idx |= (char_value << ((m - j) << 1));
	}
	++array[array_idx];
	return;
}

static R_INLINE int get_first_k_mer_index(char const * const inseq, int const k, int invalid_value)
{
	/*
	 * counts one k-mer at beginning of inseq
	 * Returns index (for later translation)
	 */

	int j, char_value, m, index;

	/*
	 * if(k>max_k)
	 * 	error("[do_count_Kmers] k must be <= %u!", max_k);
	 */

	m = k - 1;
	index = 0;
	for(j=0; j<k; ++j)
	{
		/* 	Throw error when traversing '\0' */
		if(inseq[j] == '\0')
			error("[count_Kmer] Found string terminating NULL!");

		char_value = ACGT[ (unsigned)inseq[j] ];

		if(char_value == zval)
			return invalid_value;

		/*
		 * Index reverse order (because dna Kmers do it this way)
		 * (2* ((k-1)-j))= (m-j)<<1
		 */
		index |= (char_value << ((m-j) << 1));
	}
	return index;
}

static R_INLINE void count_gc_content(char const * const inseq, int * const array, int seqlen)
{
	int i, gc;

	gc = 0;
	for(i=0; i<seqlen; ++i)
		gc += GCZ[ (unsigned)inseq[i] ];

	gc *= 100;
	/* Rprintf("[count_gc_content] gc=%u\tindex=%u\n", gc, (unsigned)(gc/seqlen)); */
	++array[ (unsigned)(gc/seqlen) ];
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		count_Kmers
 * 		k		:	Kmer length
 * 		pSeq	:	Vector of DNA sequences in Kmers are counted
 * 		pWidth	:	Vector of frame width's (number of counted frames)
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


SEXP count_Kmers(SEXP pSeq, SEXP pK, SEXP pWidth, SEXP pNn)
{
	int nRows, nCols, array_size, i;
	int k, nProtected;
	int *array;
	char *buf;

	SEXP pArray, dim, col_names, mat_dim_names;

	if(TYPEOF(pSeq) != STRSXP)
		error("[count_Kmers] pSeq must be String!");

	if(TYPEOF(pK) != INTSXP)
		error("[count_Kmers] pK must be INT!");

	if(TYPEOF(pWidth) != INTSXP)
		error("[count_Kmers] pWidth must be INT!");

	if(TYPEOF(pNn) != INTSXP)
		error("[count_Kmers] pNn must be INT!");

	nCols = LENGTH(pSeq);

	if(LENGTH(pWidth) != nCols)
		error("[count_Kmers] pSeq and pWidth must have equal length!");

	if(LENGTH(pNn) != nCols)
		error("[count_Kmers] pNn and pSeq must have equal length!");

	/* k = length of Kmer */
	k = INTEGER(pK)[0];
	if(k > max_k)
		error("[count_Kmers] k must be <= %u!", max_k);


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array which contains counted values
	 * which is organized as matrix when nCols>1
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	nRows = (1 << (k<<1));		/* 4^k */
	array_size = nRows * nCols;

	PROTECT(pArray = allocMatrix(INTSXP, nRows, nCols));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Imbue matrix structure on array
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* Column names */
	col_names = PROTECT(allocVector(STRSXP, nCols));
	buf = R_alloc(rbuf_size, sizeof(char));
    for(i=0; i<nCols; ++i)
    {
    	sprintf(buf, "%i", i + 1);
    	SET_STRING_ELT(col_names, i, mkChar(buf));
    }

	/*  Row names (DNA k-mers) */
    mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0L, create_dna_k_mers(k));
	SET_VECTOR_ELT(mat_dim_names, 1L, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);
	nProtected=4;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Do k-mer counting
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	array = INTEGER(pArray);
	memset(array, 0, ( (unsigned)array_size ) * sizeof(int) );

	for(i=0; i<nCols; ++i)
	{
		char const * const seq = CHAR(STRING_ELT(pSeq, i));
		if(!do_countCheck_Kmers(seq, array + (i * nRows), INTEGER(pNn) + i, k, (INTEGER(pWidth)[i])) )
			error("[count_Kmers] character mismatch!");
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * return
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	UNPROTECT(nProtected);
	return pArray;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		count_dna_Kmers
 * 		pSeq	:	Vector of DNA sequences in Kmers are counted
 * 		pStart	:	1-based positions of leftmost counted nucleotide
 * 		pWidth	:	Vector of frame width's (number of counted frames)
 * 		pK		:	Kmer length (One single value)
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



SEXP count_dna_Kmers(SEXP pSeq, SEXP pStart, SEXP pWidth, SEXP pK, SEXP pNn)
{
	int k, nProtected, nCols, nRows, array_size, i;
	int *array;
	SEXP pArray, dim, col_names, mat_dim_names;
	char *buf;

	/* DNA sequence		*/
	char const * const dna=CHAR(STRING_ELT(pSeq, 0));


	if(TYPEOF(pSeq) != STRSXP)
		error("[count_dna_Kmers] pSeq must be String!");

	if(LENGTH(pSeq) > 1)
		error("[count_dna_Kmers] pSeq must have length 1!");

	if(TYPEOF(pStart) != INTSXP)
		error("[count_dna_Kmers] pStart must be INT!");

	if(TYPEOF(pWidth) != INTSXP)
		error("[count_dna_Kmers] pWidth must be INT!");

	nCols = LENGTH(pStart);
	if(LENGTH(pWidth) != nCols)
		error("[count_dna_Kmers] pStart and pWidth must have equal length!");

	if(TYPEOF(pNn) != INTSXP)
		error("[count_dna_Kmers] pNn must be INT!");

	if(LENGTH(pNn) != nCols)
		error("[count_dna_Kmers] pNn and pStart must have equal length!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * k-value (number of nucleotides in DNA-motif)
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pK) != INTSXP)
		error("[count_dna_Kmers] pK must be INT!");

	if(LENGTH(pK) > 1)
		error("[count_dna_Kmers] pK must have length 1");

	k = INTEGER(pK)[0];
	if(k < 1)
		error("[count_dna_Kmers] k must be positive!");

	if(k > max_k)
		error("[count_dna_Kmers] k must be <= %u!", max_k);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array which contains counted values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nRows = (1<< (k << 1) );		/*  4^k */
	array_size = nRows * nCols;
	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attibute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column names
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(col_names = allocVector(STRSXP, nCols));
	buf = R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nCols; ++i)
	{
		sprintf(buf, "%i", INTEGER(pStart)[i]);
		SET_STRING_ELT(col_names, i, mkChar(buf));
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names: DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, create_dna_k_mers(k));
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);
	nProtected= 4 ;


	array = INTEGER(pArray);
	memset(array, 0, ((unsigned long )array_size)* sizeof(int));

	for(i=0; i<nCols; ++i)
	{
		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Add column to result data.frame
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(!do_countCheck_Kmers(dna + INTEGER(pStart)[i] - 1,
				array + (i* nRows), INTEGER(pNn) + i, k, INTEGER(pWidth)[i]))
		{
			error("[count_dna_Kmers] character mismatch!");
		}
	}

	UNPROTECT(nProtected);
	return pArray;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		rev_count_Kmers
 * 		k		:	Kmer length
 * 		inseq	:	Pointer to DNA sequence in which Kmers are counted
 * 		array	:	Array where motifs are counted
 * 		npos	:	Number of start positions of Kmers
 *
 * 		Counts complement nucleotides in reverse order
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static R_INLINE bool do_rev_count_Kmers(char const * const inseq, int * const array, int * nn, int const k, int npos)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * iter traverses inseq from right to left for npos
	 * positions
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	char const 	*c, *iter;
	unsigned 	char_value, array_idx, val_ok;
	int 		i, j, m; /* cannot be unsigned!	*/

	iter = inseq;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * k= length of Kmer; expected size of array: 4^k
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(k >max_k)
		error("[do_rev_count_Kmers] k must be <= %u!", max_k);

	m = k - 1;
	array_idx = 0;
	val_ok = 0;
	for(i=0; i<npos; ++i)
	{
		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Calculate index
		 * c walks k positions to left starting at iter
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		c = iter;
		for(j=m; j>=0; --j)
		{
			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Rprintf("[do_rev_count_Kmers] i=%u\tj=%u\t'%s'\tc=%c\n", i, j, iter, * c);
			 * Throw error when traversing '\0' !
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			if(*c == '\0')
				error("[do_rev_count_Kmers] Found string terminating NULL!");

			char_value = revACGT[(unsigned)(*c)];
			if(char_value == zval)
			{
				/*
				 * Rprintf("[do_rev_count_Kmers] zval: c=%c\tchar_val=%u\tACGT=%u\tj=%i\n", * c, char_value, ACGT[(unsigned)(* c)], j);
				 */
				if(ACGTN[ (unsigned)(*c) ] == nval)
				{
					j -= k;
					++(*nn);
					val_ok = 0;
					break;
				}
				else
				{
					/*
					 * Rprintf("[do_countCheck_Kmers] Error : j: %u\t iter: '%s'\n", j, iter);
					 */
					return false;
				}
			}
			/* Index reverse order (because dna Kmers do it this way) */
			array_idx |= (char_value << ( j<<1 ) );
			--c;
			val_ok = 1;
		}
		if(val_ok)
			++array[array_idx];
		--iter;
		array_idx = 0;
	}
	return true;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		count_dna_Kmers
 * 		pSeq	:	One DNA sequence in which Kmers are counted
 * 		pStart	:	1-based positions of rightmost counted nucleotide
 * 		pWidth	:	Vector of frame width's (number of counted frames)
 * 		pK		:	Kmer length (One single value)
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

SEXP rev_count_dna_Kmers(SEXP pSeq, SEXP pStart, SEXP pWidth, SEXP pK, SEXP pNn)
{
	int nProtected;
	int i, k, nRows, nCols, array_size;
	char *buf;
	int *array;

	SEXP pArray, dim, col_names, mat_dim_names;
	char const *const dna = CHAR(STRING_ELT(pSeq, 0));

	if(TYPEOF(pSeq) != STRSXP)
		error("[rev_count_dna_Kmers] pSeq must be String!");

	if(LENGTH(pSeq) > 1)
		error("[rev_count_dna_Kmers] pSeq must have length 1!");

	if(TYPEOF(pStart) != INTSXP)
		error("[rev_count_dna_Kmers] pStart must be INT!");

	if(TYPEOF(pWidth) != INTSXP)
		error("[rev_count_dna_Kmers] pWidth must be INT!");

	if(TYPEOF(pK) != INTSXP)
		error("[rev_count_dna_Kmers] pK must be INT!");

	if(TYPEOF(pNn) != INTSXP)
		error("[rev_count_dna_Kmers] pNn must be INT!");

	nCols = LENGTH(pStart);

	if(LENGTH(pWidth) != nCols)
		error("[rev_count_dna_Kmers] pStart and pWidth must have equal length!");

	if(LENGTH(pNn) != nCols)
		error("[rev_count_dna_Kmers] pNn and pStart must have equal length!");

	k = (int) INTEGER(pK)[0];
	if(k < 1)
		error("[rev_count_dna_Kmers] k must be positive!");

	/* k = length of Kmer */
	if(k > max_k)
		error("[rev_count_dna_Kmers] k must be <= %u!", max_k);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array which contains counted values
	 * which is organized as matrix
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	nRows=(1 << (k<<1));					/* 4^seqlen	*/
	array_size = nRows * nCols;				/* nStart	*/
	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
	array = INTEGER(pArray);
	memset(array, 0, ((unsigned long int) array_size) * sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(dim = allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column names
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(col_names = allocVector(STRSXP, nCols));
	buf=R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nCols; ++i)
	{
		sprintf(buf, "%i", i+1);
		SET_STRING_ELT(col_names, i, mkChar(buf));
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names: DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(mat_dim_names = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, create_dna_k_mers(k));
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);
	nProtected = 4;

	/* For each sequence	*/
	for(i=0; i<nCols; ++i)
	{
		/* Smallest reading index must not be <0  */
		if(INTEGER(pWidth)[i] + ((int)k) > INTEGER(pStart)[i])
			error("[rev_count_dna_Kmers] npos + k > start (i=%u)!", i);
		if(!do_rev_count_Kmers(dna+INTEGER(pStart)[i]-1, array+(i* nRows), INTEGER(pNn)+i, k, INTEGER(pWidth)[i]))

			error("[rev_count_dna_Kmers] character mismatch!");
	}
	UNPROTECT(nProtected);
	return pArray;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		count_genome_Kmers
 * 		pSeq	:   Vector of pointers to all used DNA sequences
 * 		pSeqid	:	1-based index to DNA sequence for each counting frame
 * 		pLstart	:	1-based position of leftmost counted nucleotide
 * 		pWidth	:	Number of nucleotides to be counted (frame width)
 * 		pStrand	:	1 for (+) strand, any other for (-) strand
 * 		k		:	Kmer length
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



SEXP count_genome_Kmers(SEXP pSeq, SEXP pSeqid, SEXP pLstart, SEXP pWidth, SEXP pStrand, SEXP pK, SEXP pNn)
{
	int i, k, nRows, nCols, nProtected;
	int *array;
	char *buf;

	SEXP pArray, dim, col_names, mat_dim_names;

	if(TYPEOF(pSeq) != STRSXP)
		error("[count_genome_Kmers] pSeq must be String!");

	if(TYPEOF(pLstart) != INTSXP)
		error("[count_genome_Kmers] pLstart must be INT!");

	if(TYPEOF(pWidth) != INTSXP)
		error("[count_genome_Kmers] pWidth must be INT!");

	if(TYPEOF(pStrand) != INTSXP)
		error("[count_genome_Kmers] pStrand must be INT!");

	if(TYPEOF(pK) != INTSXP)
		error("[count_genome_Kmers] pK must be INT!");

	if(TYPEOF(pNn) != INTSXP)
		error("[count_genome_Kmers] pNn must be INT!");

	nCols = LENGTH(pLstart);
	if(LENGTH(pWidth) != nCols)
		error("[count_genome_Kmers] pLstart and pWidth must have equal length!");

	if(LENGTH(pStrand) != nCols)
		error("[count_genome_Kmers] pLstart and pStrand must have equal length!");

	k = INTEGER(pK)[0];
	if(k < 1)
		error("[count_genome_Kmers] k must be positive!");

	/* k = length of Kmer */
	if(k > (int)max_k)
		error("[count_genome_Kmers] k must be <= %u!", max_k);

	if(LENGTH(pNn) != nCols)
		error("[count_genome_Kmers] pNn and pLstart must have equal length!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array which contains counted values
	 * (matrix)
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nRows= (1 << (k << 1));				/* 4^seqlen	*/
	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
	array = INTEGER(pArray);
	memset(array, 0, ((size_t)nRows) * ((size_t)nCols) * sizeof(int));

	if(nCols > 1)
	{
		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Dim attribute
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		PROTECT(dim = allocVector(INTSXP, 2));
		INTEGER(dim)[0] = nRows;
		INTEGER(dim)[1] = nCols;
		setAttrib(pArray, R_DimSymbol, dim);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Column names
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		PROTECT(col_names = allocVector(STRSXP, nCols));
		buf = R_alloc(rbuf_size, sizeof(char));
		for(i=0; i<nCols; ++i)
		{
			sprintf(buf, "%i", i + 1);
			SET_STRING_ELT(col_names, i, mkChar(buf));
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Row names: DNA k-mers
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		PROTECT(mat_dim_names = allocVector(VECSXP, 2));
		SET_VECTOR_ELT(mat_dim_names, 0L, create_dna_k_mers(k));
		SET_VECTOR_ELT(mat_dim_names, 1L, col_names);
		setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);
		nProtected = 4;
	}
	else
	{
		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Solely add DNA k-mers as names
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		setAttrib(pArray, R_NamesSymbol, create_dna_k_mers(k));
		nProtected = 1;
	}


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * For each sequence
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	for(i=0; i<nCols; ++i)
	{
		/* pSeqid index should be used as 1-based value	*/
		char const * const dna = CHAR(STRING_ELT(pSeq, INTEGER(pSeqid)[i] - 1));

		/* Do count Kmers	*/
		if(INTEGER(pStrand)[i] == 1) 	/* (+) strand */
		{
			if(!do_countCheck_Kmers(dna+INTEGER(pLstart)[i] - 1, array + (i * nRows),
					INTEGER(pNn) + i, k, INTEGER(pWidth)[i]))
			{
				error("[count_genome_Kmers] character mismatch at position %u!", i);
			}
		}
		else
		{
			if(!do_rev_count_Kmers(dna+INTEGER(pLstart)[i] + INTEGER(pWidth)[i] - 2,
					array + (i * nRows), INTEGER(pNn) + i, k, INTEGER(pWidth)[i]))
			{
				error("[count_genome_Kmers] character mismatch at position %u!");
			}
		}
	}

	UNPROTECT(nProtected);
	return pArray;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		count_splice_Kmers
 * 		pSeq	:   Vector of pointers to all used DNA sequences
 * 		pSeqid	:	1-based index to DNA sequence for each counting frame
 * 		pLend	:	1-based position of (left)  last  exon nucleotide
 * 		pRstart	:	1-based position of (right) first exon nucleotide
 * 		pWidth	:	Number of nucleotides to be counted (frame width)
 * 		pStrand	:	1 for (+) strand, any other for (-) strand
 * 		k		:	Kmer length
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

SEXP count_splice_Kmers(SEXP pSeq, SEXP pSeqid, SEXP pLend, SEXP pRstart, SEXP pWidth, SEXP pStrand, SEXP pK, SEXP pNn)
{
	int nRows, nCols, i, nProtected, nPlusMm, nMinusMm;
	size_t array_size;
	int k;
	int * array;
	char * buf;

	SEXP pArray, dim, col_names, mat_dim_names;

	if(TYPEOF(pSeq) != STRSXP)
		error("[count_splice_Kmers] pSeq must be String!");

	if(TYPEOF(pSeqid) != INTSXP)
		error("[count_splice_Kmers] pSeqid must be INT!");

	if(TYPEOF(pLend) != INTSXP)
		error("[count_splice_Kmers] pLstart must be INT!");

	if(TYPEOF(pRstart) != INTSXP)
		error("[count_splice_Kmers] pRstart must be INT!");

	if(TYPEOF(pWidth) != INTSXP)
		error("[count_splice_Kmers] pWidth must be INT!");

	if(TYPEOF(pStrand) != INTSXP)
		error("[count_splice_Kmers] pStrand must be INT!");

	if(TYPEOF(pK) != INTSXP)
		error("[count_splice_Kmers] pK must be INT!");

	if(TYPEOF(pNn) != INTSXP)
		error("[count_splice_Kmers] pNn must be INT!");


	nCols = LENGTH(pLend);
	if(LENGTH(pRstart) != nCols)
		error("[count_splice_Kmers] pLend and pRstart must have equal length");

	if(LENGTH(pWidth) != nCols)
		error("[count_splice_Kmers] pLend and pWidth must have equal length!");

	if(LENGTH(pStrand) != nCols)
		error("[count_splice_Kmers] pLstart and pStrand must have equal length!");

	if(LENGTH(pNn) != nCols)
		error("[count_splice_Kmers] pNn and pLend must have equal length!");

	k=INTEGER(pK)[0];
	if(k < 1)
		error("[count_splice_Kmers] k must be positive!");

	/* k = length of Kmer */
	if(k > max_k)
		error("[count_splice_Kmers] k must be <= %u!", max_k);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array which contains counted values
	 * (matrix)
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nRows=(1<<(k<<1));						/* 4^seqlen */
	array_size=(size_t) (nRows* nCols);		/* nStart	*/
	pArray=PROTECT(allocMatrix(INTSXP, nRows, nCols));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(dim=allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column names
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(col_names = allocVector(STRSXP, nCols));
	buf=R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nCols; ++i)
	{
		sprintf(buf, "%i", i + 1);
		SET_STRING_ELT(col_names, i, mkChar(buf));
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names: DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(mat_dim_names=allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, create_dna_k_mers(k));
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);
	nProtected = 4;

	array=INTEGER(pArray);
	memset(array, 0, array_size* sizeof(int));
	nPlusMm = 0;
	nMinusMm = 0;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * For each sequence:
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	for(i=0; i<nCols; ++i)
	{
		/*  pSeqid index should be used as 1-based value */
		char const * const dna = CHAR(STRING_ELT(pSeq, INTEGER(pSeqid)[i]-1));
		if(INTEGER(pStrand)[i] == 1) /* (+) strand	*/
		{
			/*
			 * 	pLend: 1-based, char * inseq: 0-based
			 *
			 * Window position: Rightmost nuc of rightmost k-mer inside window = rightmost nuc of left exon
			 *
			 * 			A A A A C C C C G T C C C C A G C C C C
			 * 			    x x x x									k=3, width=4
			 *
			 *
			 */
			if(!do_countCheck_Kmers(dna + INTEGER(pLend)[i] -
						INTEGER(pWidth)[i] - k + 1, array +(i * nRows), INTEGER(pNn) + i, k, INTEGER(pWidth)[i]))
			{ ++nPlusMm; }
		}
		else
		{
			if(!do_rev_count_Kmers(dna + INTEGER(pRstart)[i] +
					INTEGER(pWidth)[i] - 2, array + (i * nRows), INTEGER(pNn) + i, k, INTEGER(pWidth)[i]))
			{ ++nMinusMm; }
		}
	}

	Rprintf("[count_splice_Kmers] Finished. char mismatches: %u on (+) and %u on (-) strand.\n", nPlusMm, nMinusMm);
	UNPROTECT(nProtected);
	return pArray;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  	count_fastq_Kmers
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


fqParser * r_do_init_faq(const char* filename, int buf_capacity)
{
	fqParser * fqp=fqp_init(filename, buf_capacity);
	if(!fqp)
	{
		Rprintf("\n[fastqq] fqp_init returned 0!\n");
		return 0;
	}

	if(!fqpIsOpen(fqp))
	{
		Rprintf("\n[fastqq] Can not open file '%s'!\n", filename);
		fqp_destroy(fqp);
		return 0;
	}
	fqp_fill_rbuf(fqp);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Returns fully initalized object.
	 * - File open and not empty
	 * - rfc array countains initial filling.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	return fqp;

}

SEXP getIndexVector(long int n, int offset)
{
	long int i, off;
	char * buf;
	SEXP res;

	res=PROTECT(allocVector(STRSXP, n));
	buf=R_alloc(rbuf_size, sizeof(char));
	off=(long int) offset;
	for(i=0;i<n;++i)
	{
		sprintf(buf, "%li", i+off);
		SET_STRING_ELT(res, i, mkChar(buf));
	}
	UNPROTECT(1);
	return res;
}


SEXP createKmerCountArray(int k, int nCols, int * nRows)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create matrix count values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	size_t array_size;
	SEXP pArray, dim, col_names, mat_dim_names;

	*nRows = (1<< (k<<1) );			/* 4^k */
	array_size = (size_t) ((*nRows)* nCols);
	pArray = PROTECT(allocMatrix(INTSXP, *nRows, nCols));
	memset(INTEGER(pArray), 0, array_size * sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim=PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = *nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim names
	 * Col names	: 1:nCols
	 * Row names	: DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	PROTECT(col_names=getIndexVector(nCols, 1));
	mat_dim_names=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, PROTECT(create_dna_k_mers(k)));
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	UNPROTECT(5);
	return pArray;
}

SEXP getNucCountArray(int nCols)
{
	int i, nRows;
	size_t array_size;
	int *array;
	SEXP pArray, dim, col_names, row_names, mat_dim_names;
	char *buf;

	if(nCols<2)
		error("[getNucCountArray] nCols must be >=2!");

	nRows = nIupac;
	array_size = (size_t)(nRows * nCols);
	PROTECT(pArray = allocMatrix(INTSXP, nRows, nCols));
	array = INTEGER(pArray);
	memset(array, 0, array_size * sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim names
	 * Col names	: 1:nCols
	 * Row names	: IUPAC characters
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	col_names=PROTECT(getIndexVector(nCols, 1));
	row_names=PROTECT(allocVector(STRSXP, nRows));
	buf=R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nRows; ++i)
	{
		sprintf(buf, "%c", iupac_chars[i]);
		SET_STRING_ELT(row_names, i, mkChar(buf));
	}
	PROTECT(mat_dim_names=allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, row_names);
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	UNPROTECT(5);
	return pArray;
}

SEXP getSeqLenCountArray(int maxFastqSeqLen, int nFastqFiles)
{
	int nRows, nCols;
	size_t array_size;
	int * array;
	SEXP pArray, dim, row_names, col_names, mat_dim_names;

	nRows = maxFastqSeqLen;
	nCols = nFastqFiles;
	array_size = (size_t) (nRows * nCols);

	pArray=PROTECT(allocMatrix(INTSXP, nRows, nCols));
	array=INTEGER(pArray);
	memset(array, 0, array_size * sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim=PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim names
	 * Col names	: 1:nCols
	 * Row names	: 1:nRows
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	row_names = PROTECT(getIndexVector(nRows, 1));
	col_names = PROTECT(getIndexVector(nCols, 1));
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, row_names);
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	UNPROTECT(5);
	return pArray;
}

SEXP getSeqQualCountArray(int nCols)
{
	int nRows;
	size_t array_size;
	int *array;
	SEXP pArray, dim, col_names, row_names, mat_dim_names;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Phred counting starts at 0: 0:93		-> 94
	 * nCols = maxSeqLen
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nRows = maxPhred + 1;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	array_size = (size_t) (nRows * nCols);
	pArray=PROTECT(allocMatrix(INTSXP, nRows, nCols));
	array=INTEGER(pArray);
	memset(array, 0, array_size * sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim=PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim names:
	 * Col names	: 1:nCols
	 * Row names	: 1:nRows
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	col_names = PROTECT(getIndexVector(nCols, 1));
	row_names = PROTECT(getIndexVector(nRows, 0));
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, row_names);
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	UNPROTECT(5);
	return pArray;
}

SEXP getSeqLenArray(int nFiles)
{
	int nRows, nCols;
	size_t i, array_size;
	SEXP pArray, dim, row_names, col_names, mat_dim_names;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Array for MinSeqLen and Max SeqLen
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nRows = 2;			/* = Number of parameters	*/
	nCols = nFiles;
	array_size = (size_t) (nRows * nCols);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array and fill with default
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
	for(i=0; i<array_size; ++i)
	{
		if(i%2 == 0)
			INTEGER(pArray)[i] = 0x40000000;  /* =2^30 = 1.073.741.824 */
		else
			INTEGER(pArray)[i] = 0;			/* maxSeqLen (?)		*/
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names: minSeqLen + maxSeqLen
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	row_names = PROTECT(allocVector(STRSXP, nRows));
	SET_STRING_ELT(row_names, 0, mkChar("minSeqLen"));
	SET_STRING_ELT(row_names, 1, mkChar("maxSeqLen"));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column names	: 1:nCols
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	col_names = PROTECT(getIndexVector(nCols, 1));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Set dim names
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, row_names);
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	UNPROTECT(5);
	return pArray;
}

SEXP getGCcontentArray(int gcc_nRows, int gcc_nCols)
{
	unsigned long gcc_array_size;
	SEXP pGccRowNames, pGccColNames, pGCcount;
	SEXP dim, mat_dim_names;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare counts of GC content
	 * nRows	= maxSeqLen
	 * nCols	= nFiles
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	gcc_array_size = (unsigned long) (gcc_nRows* gcc_nCols);
	pGccRowNames = PROTECT(getIndexVector(gcc_nRows, 0));
	pGccColNames = PROTECT(getIndexVector(gcc_nCols, 1));
	pGCcount = PROTECT(allocMatrix(INTSXP, gcc_nRows, gcc_nCols));
	memset(INTEGER(pGCcount), 0, gcc_array_size* sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = gcc_nRows;
	INTEGER(dim)[1] = gcc_nCols;
	setAttrib(pGCcount, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim names
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, pGccRowNames);
	SET_VECTOR_ELT(mat_dim_names, 1, pGccColNames);
	setAttrib(pGCcount, R_DimNamesSymbol, mat_dim_names);

	UNPROTECT(5);
	return pGCcount;
}


SEXP count_fastq(SEXP pInfile, SEXP pK)
{
	int i, j, k, maxSeqLen, nFiles;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Kmer count arrays
	 * Values will be set in function call
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	int		nKmerCols;
	int		nKmerRows;
	SEXP	pKmerArray, pFirstKmerArray;
	int		* kmer_array, * first_kmer_array;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * seqLenCount
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP pSeqLenCount;
	int * seqlen_array;
	int seqlen_nRows;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Min and max seqlen
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP pSeqLen;
	int * mSeqLen_array;
	int mSeqLen_nRows;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Read count
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP pReadCount;
	int * readCount;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * N counts
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP pNCount;
	int * nNCount;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * GC content counts (Static size : 0:100 !!)
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	static const int gcc_nRows = 101;
	SEXP pGCcount;
	int * gcCount;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Nucleic acids counts
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP nacList;
	int * nucCount;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Phred counts
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP phredList;
	int * phredCount;
	const int nPhredRows = maxPhred+1;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Fastq parser
	 * seq= pointer for DNA sequence and phred qualities
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	fqParser * fqp;
	int seqlen,phred_value,nResize;
	char * seq;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Return values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP pMaxSeqLen,pnFiles,res;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Check incoming arguments
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pInfile) != STRSXP)
		error("[fastqq] pInfile must be a string!");
	nFiles = LENGTH(pInfile);

	if(TYPEOF(pK) != INTSXP)
		error("[fastqq] pK must be INT!");
	k = INTEGER(pK)[0];
	if(k<1)
		error("[fastqq] k must be positive!");

	/* k = length of Kmer */
	if(k>max_k)
		error("[fastqq] k must be <= %u!",max_k);

	maxSeqLen = default_fastq_max_seqlen;
	if(maxSeqLen<k)
		error("[fastqq] Maximal Sequence must be >=%i!",k);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare Kmer and firstKmer count arrays
	 * nRows	= 4^k
	 * nCols	= nFiles
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nKmerCols = nFiles;
	pKmerArray = PROTECT(createKmerCountArray(k, nKmerCols, &nKmerRows));
	kmer_array = INTEGER(pKmerArray);
	pFirstKmerArray = PROTECT(createKmerCountArray(k, nKmerCols, &nKmerRows));
	first_kmer_array = INTEGER(pFirstKmerArray);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare seqlen counting
	 * nRows	= maxSeqLen
	 * nCols	= nFiles
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pSeqLenCount = PROTECT(getSeqLenCountArray(maxSeqLen, nFiles));
	seqlen_nRows = maxSeqLen;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare array for MinSeqLen and MaxSeqLen
	 * nRows	= 2
	 * nCols	= nFiles
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pSeqLen = PROTECT(getSeqLenArray(nFiles));
	mSeqLen_nRows = 2;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare read counts
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pReadCount = PROTECT(allocVector(INTSXP, nFiles));
	readCount = INTEGER(pReadCount);
	memset(readCount, 0, ((size_t)nFiles)* sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare N counts
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pNCount = PROTECT(allocVector(INTSXP, nFiles));
	nNCount = INTEGER(pNCount);
	memset(nNCount, 0, ((size_t)nFiles)* sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare coutns for GC content
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pGCcount = PROTECT(getGCcontentArray(gcc_nRows, nFiles));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare Nucleic acid counts: nac
	 * nRows	= nIupac
	 * nCols	= maxSeqLen
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nacList = PROTECT(allocVector(VECSXP, nFiles));
	for(i=0;i<nFiles;++i)
		SET_VECTOR_ELT(nacList, i, PROTECT(getNucCountArray(maxSeqLen)));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare Position-wise counting of Phred quality
	 * nRows	= maxPhred + 1
	 * nCols	= maxSeqLen
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	phredList = PROTECT(allocVector(VECSXP, nFiles));
	for(i=0;i<nFiles;++i)
		SET_VECTOR_ELT(phredList, i, PROTECT(getSeqQualCountArray(maxSeqLen)));


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Open file and read sequence
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nResize = 0;

	for(i=0;i<nFiles;++i)
	{
		Rprintf("[fastqq] File (%2u/%u) '%s'", i+1, nFiles, CHAR(STRING_ELT(pInfile, i)));
		kmer_array = INTEGER(pKmerArray)+(i* nKmerRows);
		seqlen_array = INTEGER(pSeqLenCount)+(i* seqlen_nRows);
		mSeqLen_array = INTEGER(pSeqLen)+(i* mSeqLen_nRows);
		gcCount = INTEGER(pGCcount)+(i* gcc_nRows);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Initialize fastq parser structure
		 * (open file, init buffers)
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		fqp = r_do_init_faq(CHAR(STRING_ELT(pInfile, i)), default_fq_buf_capacity);

		if(!fqp)
			break;

		while(!(fqpEmpty(fqp) || fqpFastqError(fqp)))
		{
			fqp_procNuc(fqp);
			if(fqpFastqError(fqp))
			{
				Rprintf("\n[fastqq] Fastq format error:\tFile %u Read %u!\n", i+1, fqp->nSeq);
				break;
			}

			seqlen = fqp->seqlen;
			seq = fqp->pbuf;

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Resize arrays when maxSeqLen is exceeded
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			if(seqlen>maxSeqLen)
			{
				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * Resize all arrays and reset pointers
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				maxSeqLen = seqlen> 2* maxSeqLen ? seqlen : 2* maxSeqLen;
				++nResize;

				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * seqlen count
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				seqlen_nRows = maxSeqLen;
				pSeqLenCount = PROTECT(enlarge_int_mat(pSeqLenCount, seqlen_nRows, nFiles));
				seqlen_array = INTEGER(pSeqLenCount)+(i* seqlen_nRows);

				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * Nucleic acid counts
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				for(j=0;j<nFiles;++j)
					SET_VECTOR_ELT(nacList, j, PROTECT(enlarge_int_mat(VECTOR_ELT(nacList, j), nIupac, maxSeqLen)));

				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * Phred count
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				for(j=0;j<nFiles;++j)
					SET_VECTOR_ELT(phredList, j, PROTECT(enlarge_int_mat(VECTOR_ELT(phredList, j), maxPhred+1, maxSeqLen)));
			}

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Check minSeqLen and MaxSeqLen
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			if(seqlen<mSeqLen_array[0])
				mSeqLen_array[0] = seqlen;
			if(seqlen>mSeqLen_array[1])
				mSeqLen_array[1] = seqlen;

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Count sequece length
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			if(seqlen<=maxSeqLen)
				++seqlen_array[seqlen-1];

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Do kmer counting (No errors when non-nucs are found)
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */

			if(do_count_Kmers(seq, kmer_array, &(fqp->nN), k, seqlen-k+1)<0)
			{
				Rprintf("\n");
				Rprintf("[fastqq] Read sequence: '%s'\n", seq);
				error("\n[fastqq] Found string terminating 0 while Kmer count in Read %u!", fqp->nSeq);
			}
			count_first_Kmer(seq, first_kmer_array, k);

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Position dependent counting of Nucleotides
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */

			nucCount = INTEGER(VECTOR_ELT(nacList, i));
			for(j=0;j<seqlen;++j)
				++nucCount[nIupac* j+IUPAC[(unsigned char) seq[j]]];

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Count GC content
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			count_gc_content(seq, gcCount, seqlen);

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Position wise counting of phred scores
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			fqp_procPhred(fqp);
			if(fqpFastqError(fqp))
			{
				Rprintf("\n[fastqq] Fastq format error:\tFile %u Read %u!\n", i+1, fqp->nSeq);
				break;
			}
			phredCount = INTEGER(VECTOR_ELT(phredList, i));
			for(j=0;j<seqlen;++j)
			{
				/* Do regular counting */
				phred_value = (((unsigned char)seq[j])-33);
				++phredCount[(nPhredRows* j)+phred_value];
			}
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Clean up
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		INTEGER(pNCount)[i]+=fqp->nN;
		INTEGER(pReadCount)[i] = fqp->nSeq;
		if(!fqpFastqError(fqp))
			Rprintf("\tdone.\n");
		else
			Rprintf("[fastqq] File closed.\n");
		fqp_destroy(fqp);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * -> Next file
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Eventually shrink array size to actual
	 * maximal sequence length
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	int totalMaxSeqLen = 0;
	for(i=0;i<nFiles;++i)
	{
		mSeqLen_array = INTEGER(pSeqLen)+(i* mSeqLen_nRows);
		if(mSeqLen_array[1]>totalMaxSeqLen)
			totalMaxSeqLen = mSeqLen_array[1];
	}
	if(totalMaxSeqLen<maxSeqLen)
	{
		maxSeqLen = totalMaxSeqLen;
		++nResize;
		pSeqLenCount = PROTECT(cut_down_int_mat(pSeqLenCount, maxSeqLen, nFiles));
		for(j=0;j<nFiles;++j)
			SET_VECTOR_ELT(nacList, j, PROTECT(cut_down_int_mat(VECTOR_ELT(nacList, j), nIupac, maxSeqLen)));
		for(j=0;j<nFiles;++j)
			SET_VECTOR_ELT(phredList, j, PROTECT(cut_down_int_mat(VECTOR_ELT(phredList, j), maxPhred+1, maxSeqLen)));
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Construct S4 object for returning
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	pMaxSeqLen = PROTECT(allocVector(INTSXP, 1));
	INTEGER(pMaxSeqLen)[0] = maxSeqLen;

	pnFiles = PROTECT(allocVector(INTSXP, 1));
	INTEGER(pnFiles)[0] = nFiles;

	res = PROTECT(NEW_OBJECT(MAKE_CLASS("Fastqq")));
	res = SET_SLOT(res, mkChar("k"), pK);
	res = SET_SLOT(res, mkChar("maxSeqLen"), pMaxSeqLen);
	res = SET_SLOT(res, mkChar("nFiles"), pnFiles);
	res = SET_SLOT(res, mkChar("filenames"), pInfile);
	res = SET_SLOT(res, mkChar("nReads"), pReadCount);
	res = SET_SLOT(res, mkChar("nN"), pNCount);
	res = SET_SLOT(res, mkChar("seqLenCount"), pSeqLenCount);
	res = SET_SLOT(res, mkChar("gcContent"), pGCcount);
	res = SET_SLOT(res, mkChar("nac"), nacList);
	res = SET_SLOT(res, mkChar("phred"), phredList);
	res = SET_SLOT(res, mkChar("seqLen"), pSeqLen);
	res = SET_SLOT(res, mkChar("kmer"), pKmerArray);
	res = SET_SLOT(res, mkChar("firstKmer"), pFirstKmerArray);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Routine termiation
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	UNPROTECT(12+2* nFiles+nResize+2* nFiles* nResize);
	return res;
}

static R_INLINE int write_fastq_read(gzFile file, char* header, char * seq, char* phred, int seqlen, int * nWritten)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * nWitten is increased inside here,
	 * because it is also used to determine that this is
	 * not the first writing procedure
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	if((* nWritten)>0)
		gzputc(file, (int)'\n');

	if(header == 0)
	{
		// 'Empty' header given
		if(gzputs(file, "@\n") != 2)
		{
			Rprintf("[trim_fastq] Error while writing header line!\n");
			return -1;
		}
	}
	else
	{
		// Write FASTQ read header line
		gzputs(file, "@");
		gzputs(file, header);
		gzputs(file, "\n");
	}

	if(gzputs(file, seq) != seqlen)
	{
		Rprintf("[trim_fastq] Error while writing sequence: '%s' seqlen = %i!\n", seq, seqlen);
		return -1;
	}
	gzputc(file, (int)'\n');

	if(gzputs(file, "+\n") != 2)
	{
		Rprintf("[trim_fastq] Error while writing separator (+) line!\n");
		return -1;
	}
	if(gzputs(file, phred) != seqlen)
	{
		Rprintf("[trim_fastq] Error while writing quality: '%s'!\n", phred);
		return -1;
	}
	++(* nWritten);
	return 0;
}

SEXP trim_fastq(SEXP pInfile, SEXP pVals, SEXP pOutfile)
{
	int fixTrimLeft, fixTrimRight, qualTrimLeft, qualTrimRight;
	int qualDiscard, qualMask, qualMaskValue, minSeqLen;

	int qtl, qtr, qm, qualDiscardVal;
	int discard;

	gzFile keepfile, discardfile;
	fqParser * fqp;

	size_t ar_size, header_length;
	int seqlen, maxSeqLen;
	int nKept, nDiscarded;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * pointer for arrays
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	char * header, * seq, * phred;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Iterator pointers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Pointer to DNA sequence
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	char * sBegin, * sEnd, * sIter;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Pointer to phred values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	char * pBegin, * pEnd, * pIter;

	SEXP res;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Check incoming arguments
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pInfile) != STRSXP)
		error("[trim_fastq] pInfile must be a string!");
	if(LENGTH(pInfile) != 1)
		error("[trim_fastq] pInfile must have length 1!");

	if(TYPEOF(pOutfile) != STRSXP)
		error("[trim_fastq] pOutfile must be a string!");
	if(LENGTH(pOutfile) != 2)
		error("[trim_fastq] pOutfile must have length 2: keep + discard!");

	if(TYPEOF(pVals) != INTSXP)
		error("[trim_fastq] pVals must be Integer!");

	if(LENGTH(pVals) != 8)
		error("[trim_fastq] pVals must have length 8!");

	fixTrimLeft 	= INTEGER(pVals)[0];
	fixTrimRight	= INTEGER(pVals)[1];
	qualTrimLeft	= INTEGER(pVals)[2];
	qualTrimRight	= INTEGER(pVals)[3];
	qualDiscard  	= INTEGER(pVals)[4];
	qualMask		= INTEGER(pVals)[5];
	qualMaskValue	= INTEGER(pVals)[6];
	minSeqLen		= INTEGER(pVals)[7];

	qtl = qualTrimLeft+33;	/* 	The ascii value		*/
	qtr = qualTrimRight+33;
	qm = qualMask+33;
	qualDiscardVal = qualDiscard+33;

	if(fixTrimLeft<0 || fixTrimRight<0)
		error("[trim_fastq]] No negative fixTrim values allowed!");
	if(qualTrimLeft<0 || qualTrimRight<0)
		error("[trim_fastq] No negative qualTrim values allowed!");
	if(qualDiscard<0 || qualDiscard >93)
		error("[trim_fastq] qualDiscard out of range!");
	if(minSeqLen<0)
		error("[trim_fastq] minSeqLen must be positive!");
	if(qualMaskValue<33 || qualMaskValue > 127)
		error("[trim_fastq] qualMaskValue out of range (33-127)!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Open ouput files
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	keepfile = gzopen(CHAR(STRING_ELT(pOutfile, 0)), "wb");
	if(!keepfile)
		error("[trim_fastq] Could not open '%s' for writing!", CHAR(STRING_ELT(pOutfile, 0)));
	discardfile = gzopen(CHAR(STRING_ELT(pOutfile, 1)), "wb");
	if(!discardfile)
		error("[trim_fastq] Could not open '%s' for writing!", CHAR(STRING_ELT(pOutfile, 1)));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Output related values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nKept = 0;
	nDiscarded = 0;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Initialize fastq-parser structure
	 * (open file, init read buffers)
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	fqp = r_do_init_faq(CHAR(STRING_ELT(pInfile, 0)), default_fq_buf_capacity);
	if(!fqp)
		return R_NilValue;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Initalize sequence and phred buffers
	 * and iterators
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	maxSeqLen = 1000;
	ar_size = ((size_t)maxSeqLen)+1;

	seq = R_alloc(ar_size, sizeof(char));
	seq[maxSeqLen] = '\0';
	phred = R_alloc(ar_size, sizeof(char));
	phred[maxSeqLen] = '\0';

	/* Iterators */
	sBegin = seq;
	pBegin = phred;


	while(!(fqpEmpty(fqp) || fqpFastqError(fqp)))
	{
		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Fill sequence and phred buffer
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		discard = 0;

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Header line is copied (for writing to output file
		 * and fqp->r_iter pointer is left unchanged.
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(* fqp->r_iter == fqp_char_sqDelim)
		{
			header_length = find_eol(fqp->r_iter);
			if(header_length > 1)
			{
				// First character (='@') is excluded here
				// because the writing procedure
				// writes '@' from itself.
				header = R_alloc(header_length, sizeof(char));
				--header_length;
				strncpy(header, fqp->r_iter + 1, header_length);
				header[header_length] = '\0';
			}
			else
			{
				header = 0;
			}
		}
		else
			header = 0;

		fqp_procNuc(fqp);
		if(fqpFastqError(fqp))
		{
			Rprintf("\n[trim_fastq] Fastq format error: Read %u!\n", fqp->nSeq);
			break;
		}

		seqlen = fqp->seqlen;

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Resize buffers when appropriate
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(seqlen > maxSeqLen)
		{
			/* resize buffers */
			maxSeqLen = 2* maxSeqLen > seqlen ? 2* maxSeqLen : seqlen;
			ar_size = ((size_t) maxSeqLen)+1;
			seq = R_alloc(ar_size, sizeof(char));
			seq[maxSeqLen] = '\0';
			phred = R_alloc(ar_size, sizeof(char));
			phred[maxSeqLen] = '\0';
		}
		/* + + + + + + + + + + + + + + + + + + + */
		memcpy(seq, fqp->pbuf, ((unsigned long)seqlen)* sizeof(char));
		sBegin = seq;
		sEnd = seq+seqlen;
		* sEnd = '\0';

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Fill phred buffer
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		fqp_procPhred(fqp);
		if(fqpFastqError(fqp))
		{
			Rprintf("\n[trim_fastq] Fastq format error: Read %u!\n", fqp->nSeq);
			break;
		}
		memcpy(phred, fqp->pbuf, ((unsigned long) seqlen)* sizeof(char));
		pBegin = phred;
		pEnd = phred+seqlen;
		* pEnd = '\0';


		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * A) Any phred < qualDiscard ?
		 * -> Discard = Read will be written to
		 * 				discard file
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		pIter = phred;
		while(pIter != pEnd)
		{
			if(((unsigned char)* pIter)<qualDiscardVal)
			{
				discard = 1;
				break;
			}
			++pIter;
		}
		if(discard)
		{
			if(write_fastq_read(discardfile, header, seq, phred, seqlen, &nDiscarded))
			{
				/*
				 * Rprintf("[trim_fastq] Fix trim discard : Rest-len =  %i < %i\n",seqlen,minSeqLen);
				 */
				gzclose(keepfile);
				gzclose(discardfile);
				error("[trim_fastq] Write error: Read %u!\n", fqp->nSeq);
			}
			continue;
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * B)
		 * Fixed trimming:
		 * Fixed number of nucs from left and right
		 * will be removed
		 *
		 * When rest seq is shorter than minSeqLen
		 * -> Discard (write to discard file)
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		sBegin += fixTrimLeft;
		pBegin += fixTrimLeft;

		sEnd -= fixTrimRight;
		pEnd -= fixTrimRight;
		*sEnd = '\0';
		*pEnd = '\0';
		seqlen = (int)(sEnd - sBegin);

		if(seqlen<minSeqLen)
		{
			/*
			 * Rprintf("[trim_fastq] Fix trim discard : Rest-len =  %i < %i\n", seqlen, minSeqLen);
			 */
			if(write_fastq_read(discardfile, header, sBegin, pBegin, seqlen, &nDiscarded))
			{
				gzclose(keepfile);
				gzclose(discardfile);
				error("[trim_fastq] Write error: Read %u!\n", fqp->nSeq);
			}
			continue;
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * C)
		 * Quality based trim
		 * Starting form both ends, all positions are removed
		 * until first phred exceeds qualTrim value
		 *
		 * When res seq is shorter than minSeqLen
		 * -> Discard (write to discard file)
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */



		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Qual based trim on prefix
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		while(((unsigned char)* pBegin)<qtl && sBegin<sEnd)
		{
			++pBegin;
			++sBegin;
		}
		seqlen = (int)(sEnd-sBegin);


		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Qual based trim on suffix
		 * Shift one to left -> Now on last phred value
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		--sEnd;
		--pEnd;

		while( ((unsigned char) *pEnd) < qtr && sBegin <= sEnd)
		{
			--sEnd;
			--pEnd;
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Shift one to right
		 * -> Now one bedind last output phred
		 * + Insert end of string marker '\0'
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		++sEnd;
		++pEnd;
		*sEnd = '\0';
		*pEnd = '\0';
		seqlen = (int)(sEnd - sBegin);

		if(seqlen < minSeqLen)
		{
			/*
			 * Rprintf("[trim_fastq] Qual trim discard: rest-len =  %i < %i!\n", seqlen, minSeqLen);
			 */
			if(write_fastq_read(discardfile, header, sBegin, pBegin, seqlen, &nDiscarded))
			{
				gzclose(keepfile);
				gzclose(discardfile);
				error("[trim_fastq] Write error: Read %u!\n", fqp->nSeq);
			}
			continue;
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * D)
		 * Quality based mask
		 * For positions iwth phred < qualMask
		 * Nucs will be overwritten with qualMaskValue
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		sIter = sBegin;
		pIter = pBegin;
		while(sIter != sEnd)
		{
			if((unsigned char)* pIter<qm)
				* sIter = (char) qualMaskValue;
			++sIter;
			++pIter;
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * D) Write to keepfile
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		/*
		 * Rprintf("[trim_fastq] keep seq: '%s' phred: '%s'\n", sBegin, pBegin);
		 */
		if(write_fastq_read(keepfile, header, sBegin, pBegin, seqlen, &nKept))
		{
			gzclose(keepfile);
			gzclose(discardfile);
			error("[trim_fastq] Write error!");
		}

	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Clean up
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(fqpFastqError(fqp))
		Rprintf("[trim_fastq] Fastq format error at file end.\n");

	fqp_destroy(fqp);

	// Append newline
	if(nKept > 0)
		gzputs(keepfile, "\n");
	if(nDiscarded > 0)
		gzputs(discardfile, "\n");

	gzclose(keepfile);
	gzclose(discardfile);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Routine termination
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	res = PROTECT(allocVector(INTSXP, 2));
	INTEGER(res)[0] = nKept;
	INTEGER(res)[1] = nDiscarded;
	UNPROTECT(1);
	return res;
}

SEXP fastq_KmerSubset_locs(SEXP pInfile, SEXP pK, SEXP pKmerIdx)
{
	int 		nFiles;
	int 		i, j, m, p;			/* running indices	*/
	int			k;					/*  length of DNA K-mer motif	*/

	int			* kmerIdx;
	int			kArraySize;
	int			other_value;

	int			nCount;
	int  		kMerIdxLen;
	int 		* translation_vector;

	fqParser 	* fqp;
	int 		seqlen, seqpos, taix, raw_index;
	char 		* seq; 				/* pointer for DNA sequence */
	int 		* array;
	const char 	* filename;

	char		* row_names_buf;
	unsigned 	array_size;
	int 		maxSeqLen, maxFileSeqLen, nRows, nCols, nResize = 0;

	SEXP pReturnList, pArray, dim, row_names, col_names;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Check incoming arguments
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pInfile) != STRSXP)
		error("[fastq_KsubLocs] pInfile must be a string!");

	nFiles = LENGTH(pInfile);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * k = Length of DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pK) != INTSXP)
		error("[fastq_KsubLocs] pK must be INT!");

	if(INTEGER(pK)[0]<1)
		error("[fastq_KsubLocs] k must be positive!");
	if(INTEGER(pK)[0]>max_k)
		error("[fastq_KsubLocs] k must be <= %u!", max_k);
	k = INTEGER(pK)[0];

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Counted k-mers are given as vector of indices
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pKmerIdx) != INTSXP)
		error("[fastq_KsubLocs] pKmerIdx must be INT!");

	kMerIdxLen = LENGTH(pKmerIdx);
	kmerIdx = INTEGER(pKmerIdx);
	kArraySize = 1<<(k<<1); 			/* 	= 4^k;	*/

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Given k-mer indices are all within valid range for k
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	for(i=0;i<kMerIdxLen;++i)
	{
		if(kmerIdx[i]<0 || kmerIdx[i] >= kArraySize)
			error("[fastq_KsubLocs] pKmerIdx [i] %i (%i)out of range!", i, kmerIdx[i]);
	}
	if(kMerIdxLen>=kArraySize)
		error("[fastq_KsubLocs] Too many indexes given!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Static sequence position arrays
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	maxSeqLen = default_fastq_max_seqlen;
	if(maxSeqLen<k)
		error("[fastq_KsubLocs] Maximal Sequence must be >= k (%u)!", k);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create translation vector
	 * Used for modified counting algorithm:
	 * A) First "original" k-mer  indices are calculated.
	 * B) Raw indices are translated -> taix (0, 1, 2, ...)
	 * C) Translated indices represent row indices in
	 * 		result array
	 *
	 * Translation vector contains the translated value (taix)
	 * at the "original" k-mer index position:
	 * taix=translation_vector[raw_index]
	 *
	 * All other positions are filled with 'kMerIdxLen'
	 * value which accumulates coutns for everything else
	 * (including non-Nucs) into last row of result array
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	nCount = kMerIdxLen+1;
	translation_vector = (int*) R_alloc((size_t)kArraySize, sizeof(int));
	for(i=0; i<kArraySize; ++i)
		translation_vector[i]=kMerIdxLen;		/* Translated index of "other" */
	for(i=0; i<kMerIdxLen; ++i)
		translation_vector[ kmerIdx[i] ] = i;		/* Translated index for given indexes */


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Find arbitrary index which will be translated
	 * into "other".
	 * Will be used when non-nucleotide characters
	 * are found.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	i = 0;
	while(translation_vector[i] != kMerIdxLen && i<kMerIdxLen)
		++i;
	if(i == kMerIdxLen)
		error("[fastq_KsubLocs] Error while fixing other value!");
	other_value = i;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare return list
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pReturnList = PROTECT(allocVector(VECSXP, nFiles));

	row_names_buf = (char* ) R_alloc(((size_t)k)+1, sizeof(char));
	row_names_buf[k] = '\0';
	nRows = nCount;
	nResize = 0;
	for(i=0; i<nFiles; ++i)
	{

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Create counting array and add to return list
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		maxSeqLen = default_fastq_max_seqlen;
		maxFileSeqLen = 0;
		nCols = maxSeqLen-k+1;
		array_size = (unsigned) (nRows * nCols);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Construct array and add to list
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
		memset(INTEGER(pArray), 0, array_size* sizeof(int));

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Dim attribute
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		dim = PROTECT(allocVector(INTSXP, 2));
		INTEGER(dim)[0] = nRows;
		INTEGER(dim)[1] = nCols;
		setAttrib(pArray, R_DimSymbol, dim);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Dim names
		 * Col names	: 1:(maxSeqLen-k+1)
		 * Row names	: DNA k-mer motifs (+ "other")
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */

		row_names = PROTECT(allocVector(STRSXP, nCount));
		for(p=0; p<kMerIdxLen; ++p) /* nCount-1 */
		{
			for(j=0, m=(2* (k-1)); j<k; ++j)
			{
				row_names_buf[j] = (char) rev_ACGT[(kmerIdx[p] >> m) & 3];
				m -= 2;
			}
			SET_STRING_ELT(row_names, p, mkChar(row_names_buf));
		}
		SET_STRING_ELT(row_names, kMerIdxLen, mkChar("other"));

		col_names = PROTECT(getIndexVector(nCols, 1));

		SEXP mat_dim_names = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(mat_dim_names, 0, row_names);
		SET_VECTOR_ELT(mat_dim_names, 1, col_names);
		setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);


		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Initialize fastq-parser structure
		 * (open file, initialize buffers)
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		Rprintf("[fastq_KsubLocs] File (%2u/%u)", i+1, nFiles);
		array = INTEGER(pArray);
		filename = CHAR(STRING_ELT(pInfile, i));

		fqp = fqp_init(filename, default_fq_buf_capacity);
		if(!fqp)
			error("\n[fastq_KsubLocs] fqp_init returned 0!\n");
		if(!fqpIsOpen(fqp))
		{
			fqp_destroy(fqp);
			error("\n[fastq_KsubLocs] Can not open file '%s'!\n", filename);
		}
		fqp_fill_rbuf(fqp);
		Rprintf(" '%s' ", filename);


		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Process file content
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		while(!(fqpEmpty(fqp) || fqpFastqError(fqp)))
		{
			fqp_procNuc(fqp);
			if(fqpFastqError(fqp))
			{
				Rprintf("\n[fastq_KsubLocs] Fastq format error:\tFile %u Read %u!\n", i+1, fqp->nSeq);
				break;
			}

			seqlen = fqp->seqlen;
			if(seqlen>maxFileSeqLen)
				maxFileSeqLen = seqlen;

			seq = fqp->pbuf;

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Resize array when needed
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			if(seqlen>maxSeqLen)
			{
				maxSeqLen =  seqlen > (2 * maxSeqLen) ? seqlen : 2 * maxSeqLen;
				nCols = maxSeqLen -k + 1;
				pArray = PROTECT(enlarge_int_mat(pArray, nRows, nCols));
				array = INTEGER(pArray);
				++nResize;
			}

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Do K-mer counting
			 * (No errors when non-nucs are found)
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */

			for(seqpos = 0; seqpos <= (seqlen - k); ++seqpos)
			{
				raw_index = get_first_k_mer_index(seq + seqpos, k, other_value);
				taix = translation_vector[raw_index];
				++array[(seqpos* nCount)+taix];
			}

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Skip phred values
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			fqp_procPhred(fqp);
			if(fqpFastqError(fqp))
			{
				Rprintf("\n[fastq_KsubLocs] Fastq format error:\tFile %u Read %u!\n", i+1, fqp->nSeq);
				break;
			}
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Add to returned list
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(maxSeqLen>maxFileSeqLen)
		{
			nCols = maxFileSeqLen-k+1;
			pArray = PROTECT(cut_down_int_mat(pArray, nRows, nCols));
			++nResize;
		}
		SET_VECTOR_ELT(pReturnList, i, pArray);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Cleanup
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(!fqpFastqError(fqp))
			Rprintf("\tdone.\n");
		else
			Rprintf("[fastq_KsubLocs] File closed.\n");
		fqp_destroy(fqp);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * -> Next file
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	}
	UNPROTECT(1+5* nFiles+nResize);
	return pReturnList;
}


SEXP get_Kmer_Index(SEXP pSequence, SEXP pK)
{
	int n, k, i;
	const char * seq;
	SEXP res;

	if(TYPEOF(pSequence) != STRSXP)
		error("[get_Kmer_index] pSequence must be string!");
	n = LENGTH(pSequence);
	if(TYPEOF(pK) != INTSXP)
		error("[get_Kmer_index] pK must be Integer!");
	k = INTEGER(pK)[0];
	res = PROTECT(allocVector(INTSXP, n));

	for(i=0;i<n;++i)
	{
		seq = CHAR(STRING_ELT(pSequence, i));
		if(strlen(seq) != ((unsigned)k))
			error("[get_Kmer_index] All input sequences must have length k!");
		INTEGER(res)[i] = get_first_k_mer_index(seq, k, -1);
	}
	UNPROTECT(1);
	return res;
}

SEXP get_kmer(SEXP pKmerIndex, SEXP pK)
{
	int nSeq, len;
	int i, j, k;
	int * index;
	char * c;
	SEXP res;

	if(TYPEOF(pKmerIndex) != INTSXP)
		error("[get_kmer] pKmerIndex must be Int!");
	if(TYPEOF(pK) != INTSXP)
		error("[get_kmer] pK must be Int!");

	nSeq = LENGTH(pKmerIndex);
	res = PROTECT(allocVector(STRSXP, nSeq));

	index = INTEGER(pKmerIndex);

	len = INTEGER(pK)[0];
	if(len<1)
		error("[get_kmer] k must be positive!");
	if(len>max_k)
		error("[get_kmer] k must be <= max_k!");


	c = (char* ) R_alloc(((size_t)len) + 1, sizeof(char));
	c[len] = '\0';

	for(i=0; i<nSeq; ++i)
	{
		for(j=0, k=(2* (len-1));j<len;++j)
		{
			c[j]=(char) rev_ACGT[(index[i] >> k) & 3];
			k-=2;
		}
		SET_STRING_ELT(res, i, mkChar(c));
	}
	UNPROTECT(1);
	return res;
}

SEXP fastq_Kmer_locs(SEXP pInfile, SEXP pK)
{
	int			nFiles;
	int			i;					/* running indices					*/
	int 		k;					/* length of DNA K-mer motif		*/
	int			nCount;

	fqParser 	* fqp;
	int			seqlen, seqpos;
	char * seq; 				/* pointer for DNA sequence			*/
	int * array;
	const char * filename;
	int array_index;

	int nResize, nRows, nCols, maxSeqLen, maxFileSeqLen;
	unsigned array_size;
	SEXP pReturnList, pArray, dim, row_names, col_names, mat_dim_names;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Check incoming arguments
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pInfile) != STRSXP)
		error("[fastq_Klocs] pInfile must be a string!");
	nFiles = LENGTH(pInfile);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * k = length of k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pK) != INTSXP)
		error("[fastq_Klocs] pK must be INT!");
	if(INTEGER(pK)[0] < 1)
		error("[fastq_Klocs] k must be positive!");
	if(INTEGER(pK)[0] > max_k)
		error("[fastq_Klocs] k must be <= %u!", max_k);
	k = INTEGER(pK)[0];
	nCount = 1<<(k<<1);				/* 	=4^k	*/


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Static sequence position arrays
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	maxSeqLen = default_fastq_max_seqlen;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare return list
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	pReturnList = PROTECT(allocVector(VECSXP, nFiles));

	nResize = 0;
	for(i=0; i<nFiles; ++i)
	{
		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Create counting array and add to return list
		 * Row names	: DNA k-mer motifs
		 * Col names	: 1:lastCountPosition = (maxSeqLen-k+1)
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		maxSeqLen = default_fastq_max_seqlen;
		maxFileSeqLen = 0;
		nCols = maxSeqLen - k + 1;

		row_names = PROTECT(create_dna_k_mers(k));
		col_names = PROTECT(getIndexVector(nCols, 1));

		nRows = length(row_names);
		nCols = length(col_names);
		array_size = (unsigned) (nRows * nCols);

		pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
		memset(INTEGER(pArray), 0, array_size* sizeof(int));

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Dim attribute
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		dim = PROTECT(allocVector(INTSXP, 2));
		INTEGER(dim)[0] = nRows;
		INTEGER(dim)[1] = nCols;
		setAttrib(pArray, R_DimSymbol,dim);

		mat_dim_names = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(mat_dim_names, 0, row_names);
		SET_VECTOR_ELT(mat_dim_names, 1, col_names);
		setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Initialize fastq-parser structure
		 * (open file, init buffers)
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		Rprintf("[fastq_Klocs] File (%2u/%u)", i+1, nFiles);
		array = INTEGER(pArray);
		filename = CHAR(STRING_ELT(pInfile, i));

		fqp = fqp_init(filename, default_fq_buf_capacity);
		if(!fqp)
			error("\n[fastq_Klocs] fqp_init returned 0!\n");
		if(!fqpIsOpen(fqp))
		{
			fqp_destroy(fqp);
			error("\n[fastq_Klocs] Can not open file '%s'!\n", filename);
		}
		fqp_fill_rbuf(fqp);
		Rprintf(" '%s' ", filename);


		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Process file content
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		while(!(fqpEmpty(fqp) || fqpFastqError(fqp)))
		{
			fqp_procNuc(fqp);
			if(fqpFastqError(fqp))
			{
				Rprintf("\n[fastq_Klocs] Fastq format error:\tFile %u Read %u!\n", i+1, fqp->nSeq);
				break;
			}
			seqlen = fqp->seqlen;
			if(seqlen>maxFileSeqLen)
				maxFileSeqLen = seqlen;
			seq = fqp->pbuf;

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Resize array when needed
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			if(seqlen>maxSeqLen)
			{
				maxSeqLen =  seqlen>(2* maxSeqLen) ? seqlen : 2* maxSeqLen;
				nCols = maxSeqLen-k+1;
				pArray = PROTECT(enlarge_int_mat(pArray, nRows, nCols));
				array = INTEGER(pArray);
				++nResize;
			}

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 *  Do Kmer counting (No errors when non-nucs are found)
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			for(seqpos=0; seqpos <= (seqlen-k); ++seqpos)
			{
				// -1 = "other" value (returned when Non Nucs appear)
				array_index = get_first_k_mer_index(seq+seqpos, k, -1);
				if(array_index >= 0)
					++array[(seqpos * nCount) + array_index];
			}

			/* * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Skip phred values
			 * * * * * * * * * * * * * * * * * * * * * * * * * * */
			fqp_procPhred(fqp);
			if(fqpFastqError(fqp))
			{
				Rprintf("\n[fastq_Klocs] Fastq format error:\tFile %u Read %u!\n", i+1, fqp->nSeq);
				break;
			}
		}

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Add to returned list
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(maxSeqLen>maxFileSeqLen)
		{
			nCols = maxFileSeqLen - k + 1;
			pArray = PROTECT(cut_down_int_mat(pArray, nRows, nCols));
			++nResize;
		}
		SET_VECTOR_ELT(pReturnList, i, pArray);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Cleanup
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
		if(!fqpFastqError(fqp))
			Rprintf("\tdone.\n");
		else
			Rprintf("[fastq_Klocs] File closed.\n");
		fqp_destroy(fqp);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * -> Next file.
		 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	}
	UNPROTECT(1+5* nFiles+nResize);
	return pReturnList;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		Count fasta.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


SEXP write_fai(SEXP pInfile, SEXP pOutfile)
{
	char * head, * buf, * buf_iter, * head_iter, * buf_end, * head_end;

	// bytes_per_line = lineLen+1
	long int lineLen, seq_len, first_base_offset, bases_per_line;
	long int total_seq_len;
	unsigned i, nSeqs, nFiles;

	gzFile fin;
	FILE * fout;

	if(TYPEOF(pInfile) != STRSXP)
		error("[write_fai] pInfile must be a string!");

	nFiles = (unsigned)LENGTH(pInfile);

	if(TYPEOF(pOutfile) != STRSXP)
		error("[write_fai] pOutfile must be a string!");

	if((unsigned)LENGTH(pOutfile) != nFiles)
		error("[write_fai] pInfile and pOutfile must have equal length!");

	setlocale(LC_ALL, "");	/* thousands sep */

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create buffer for sequence and qualities
	 * unsigned buf_size = 1024; (stat_defs.h)
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	head = (char*) R_alloc(rbuf_size, sizeof(char));
	head_end = head + buf_size - 1;
	*head_end = '\0';
	buf=(char*)  R_alloc(rbuf_size, sizeof(char));
	buf_end = buf + buf_size - 1;
	*buf_end = '\0';

	/* bytes_per_line = lineLen+1 */
	lineLen = 0;
	seq_len = 0;
	total_seq_len = 0;
	nSeqs = 0;
	bases_per_line = 0;
	first_base_offset = 0;

	for(i=0; i<nFiles; ++i)
	{
		fin = gzopen(CHAR(STRING_ELT(pInfile, i)), "rb");
		if(fin == NULL)
			error("[write_fai] Infile (%u) '%s' does not exist!", i, CHAR(STRING_ELT(pInfile, i)));
		fout = fopen(CHAR(STRING_ELT(pOutfile, i)), "w");
		if(fout == NULL)
			error("[write_fai] Outfile (%u) '%s' cannot be opened!", i, CHAR(STRING_ELT(pOutfile, i)));

		Rprintf("[write_fai] File (%2u/%u)", i+1, nFiles);

		buf[0] = '\0';
		if(!gzgets(fin, buf, buf_size))
			break;

		while(!gzeof(fin))
		{
			if(buf[0] == '>')
			{
				++nSeqs;

				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * Copy head
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				buf_iter = buf;
				head_iter = head;
				while(!(* buf_iter == '\0' || * buf_iter == '\n' || buf_iter>=buf_end))
				{
					*head_iter = *buf_iter;
					++head_iter;
					++buf_iter;
				}
				if(head_iter<head_end)
					*head_iter = '\0';
				/*
				 * Rprintf("[write_fai] Header    line: len=%2li\tseq='%s'\n", strlen(head), head);
				 */


				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * First sequence line
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				first_base_offset = gztell(fin);
				if(!gzgets(fin, buf, buf_size))
					break;

				/* Skip comment lines */
				while(!gzeof(fin) && buf[0] == ';')
				{
					first_base_offset = gztell(fin);
					if(!gzgets(fin, buf, buf_size))
						break;
				}

				lineLen = (long) strlen(buf);

				/*
				 * Count Nucleotide characters
				 * Exclude newLine from count
				 */
				if(lineLen>0)
					--lineLen;
				seq_len = lineLen;
				bases_per_line = lineLen;

				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * Count subsequent lines
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				while(!gzeof(fin) && buf[0] != '>')
				{
					if(!gzgets(fin, buf, buf_size))
						break;
					lineLen = (long) strlen(buf);
					if(lineLen>0)
						--lineLen;
					buf[lineLen] = '\0';	// cosmetic
					if(!(buf[0] != '>' || buf[0] == ';'))
						seq_len += lineLen;
				}
			}
			/* Add Nucleotides from last line in file */
			if(gzeof(fin))
				seq_len += lineLen;
			total_seq_len += seq_len;
			fprintf(fout, "%s\t%lu\t%lu\t%lu\t%lu\n", head, seq_len, first_base_offset, bases_per_line, bases_per_line+1);
		}
		gzclose(fin);
		fclose(fout);
		Rprintf("\t%3u sequences %'12lu total length.\n", nSeqs, total_seq_len);
	}
	return R_NilValue;
}



void insertSeqName(SEXP col_names, R_xlen_t i, faTraverse * fat)
{
	char *seqname;

	if(fatNewSeq(fat))
	{
		seqname = fat_getSeqName(fat);
		SET_STRING_ELT(col_names, i, mkChar(seqname));
		/* Rprintf("[count_fasta_Kmers] New seq  : '%s'\n", seqname); */
		free(seqname);
	}
}

faTraverse * r_do_init_fat(const char * filename, int k)
{
	faTraverse * fat;
	/*
	 * Returns memory initialized structure
	 * - File may be closed or empty
	 * - pos buffer still is empty
	 */

	fat = fat_init(filename, k);

	if(!fat)
	{
		Rprintf("[count_fasta_Kmers] fat_init returned 0!\n");
		return 0;
	}

	if(!fatIsOpen(fat))
	{
		Rprintf("[count_fasta_Kmers] Can not open file '%s'!\n", filename);
		fat_destroy(fat);
		return 0;
	}

	if(fatIsEof(fat))
	{
		Rprintf("[count_fasta_Kmers] Opened file '%s' is empty!", filename);
		fat_destroy(fat);
		return 0;
	}

	Rprintf("[count_fasta_Kmers] Opened file '%s'.\n", filename);

	if(fatCheckFill(fat))
	{
		Rprintf("[count_fasta_Kmers] Buffer initialization failed.\n");
		fat_destroy(fat);
		return 0;
	}

	/*
	 * Returns a full initialized faTraverse structure:
	 * - File is opened and non-empty.
	 * - faTraverse buffer pos-array contains initial filling.
	 */
	return fat;
}


SEXP count_fasta_Kmers(SEXP pFasta, SEXP pK)
{
	int k, nCols;
	int nProtected;
	int nRows, array_size;
	char * buf;
	int i, nn = 0;
	int * array;

	int iCol;					/* column index: For inserting column names		*/
	unsigned array_offset;		/* For do_countCheck_Kmers: = (iCol-1)* nRows	*/

	faTraverse * fat;

	SEXP pArray, col_names, dim, mat_dim_names;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Check incoming args
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nProtected = 0;
	if(TYPEOF(pFasta) !=  STRSXP)
		error("[count_fasta_Kmers] pFasta must be a string!");
	if(TYPEOF(pK) != INTSXP)
		error("[count_fasta_Kmers] pK must be INT!");

	k = INTEGER(pK)[0];
	if(k < 1)
		error("[count_fasta_Kmers] k must be positive!");
	if(k > max_k)
		error("[count_fasta_Kmers] k must be <= %u!", max_k);

	/* See seqTools.h for constant */
	nCols = (int) default_fasta_Kmers_col_number;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create array which contains counted values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	nRows = (1 << (k<<1) );				/* 4^seqlen */
	array_size = (1 << (k<<1) ) * nCols;	/* nStart   */
	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
	++nProtected;
	col_names = R_NilValue;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Dim attribute
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	dim = PROTECT(allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column names:
	 * First insert numbers
	 * Then write fasta-seq id's.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	PROTECT(col_names = allocVector(STRSXP, nCols));
	++nProtected;
	buf = R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nCols; ++i)
	{
		sprintf(buf, "%i", i+1);
		SET_STRING_ELT(col_names, i, mkChar(buf));
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names: DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	++nProtected;
	SET_VECTOR_ELT(mat_dim_names, 0, create_dna_k_mers(k));
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare array and writing of indices
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	array = INTEGER(pArray);
	memset(array, 0, ((unsigned long) array_size) * sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Open file and read sequence
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	fat = r_do_init_fat(CHAR(STRING_ELT(pFasta, 0)), k);
	if(!fat)
	{
		UNPROTECT(nProtected);
		return R_NilValue;
	}
	iCol = 0;
	array_offset = 0;
	if(fat_checkNewSeq(fat))
	{
		insertSeqName(col_names, iCol, fat);
		++iCol;
		fat_skipSeqHeader(fat);
	}
	//Rprintf("[count_fasta_Kmers] Post initial fat_checkNewSeq das->r_iter: '%s'\n", fat->das->r_iter);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Traverse file and count DNA k-mers
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	while(!fatEmpty(fat))
	{
		if(fat_findKarray(fat))
		{
			//Rprintf("[count_fasta_Kmers] findKarray->fat_ok.\n");
			if(fatNucReady(fat))
			{
				/* * * * * * * * * * * * * * * * * * * * * * * * * * *
				 * Number of K-mer starting points
				 * npos=fat->das->npPos-k+1;
				 * * * * * * * * * * * * * * * * * * * * * * * * * * */
				//Rprintf("[count_fasta_Kmers] fatNucReady. array_offset: %u\tk=%i\tnpos: %i\npos:'%s'\n", array_offset, k, fat->das->npPos, fat->das->pos);
				if(!do_countCheck_Kmers(fat->das->pos, array + array_offset, &nn, k, fat->das->npPos - k + 1))
					error("[count_fasta_Kmers] character mismatch!");
			}
			if(fat_checkNewSeq(fat))
			{
				if(iCol >= nCols)
				{
					/* Rprintf("[count_fasta_Kmers] Enlarge array.\n"); */
					nCols *= 2;
					pArray = PROTECT(enlarge_int_mat(pArray, nRows, nCols));
					mat_dim_names = getAttrib(pArray, R_DimNamesSymbol);
					col_names = VECTOR_ELT(mat_dim_names, 1);

					++nProtected;
					array = INTEGER(pArray);
				}
				insertSeqName(col_names, iCol, fat);
				++iCol;
				fat_skipSeqHeader(fat);
				array_offset += (unsigned)nRows;
			}
		}
	}
	fat_destroy(fat);
	if(nn > 0)
		Rprintf("[count_fasta_Kmers] Info: Found %i N's.\n", nn);

	Rprintf("[count_fasta_Kmers] done.\n");

	pArray=PROTECT(cut_down_int_mat(pArray, nRows, iCol));
	++nProtected;
	UNPROTECT(nProtected);
	return pArray;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 	get_column_quantiles
 * 	Expects a vector of quantiles (pQuant) and
 * 	data.frame (pDf) where each column contains relative quantities (sums up to 1)
 * 	Steps down each column and writes values for each quantile
 * 	into output data.frame
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


SEXP get_column_quantiles(SEXP pQuant, SEXP pDf)
{
	unsigned nwRows, nrRows, nCols, i, j, k;
	int nProtected;
	double val, *q;
	char *buf;
	SEXP rQvec, wQvec, dflist, col_names, in_col_names, row_names;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * pQuant is expected to contain unique
	 * and ascending sorted values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pQuant) != REALSXP)
		error("[get_col_quantiles] pQuant must be REAL!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * pDf is expected to be a data.frame
	 * where all columns contain relative values
	 * (between 0 and 1, so direct comparison to pQuant
	 * works) and all columns sum up to 1.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */
	if(TYPEOF(pDf) != VECSXP)
		error("[get_col_quantiles] pDf must be VECSXP!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Nr of Rows for writing = Nr of quantile values
	 * * * * * * * * * * * * * * * * * * * * * * * * * * */

	nwRows = (unsigned) LENGTH(pQuant);

	// Nr of Rows for reading = Nr of rows in pDf
	nrRows = (unsigned) LENGTH(VECTOR_ELT(pDf, 0));

	// Nr of Columns: equal in both data.frames
	nCols = (unsigned) LENGTH(pDf);

	q = REAL(pQuant);
	nProtected = 0;

	// read and write qual-vectors
	PROTECT(dflist = allocVector(VECSXP, nCols));
	++nProtected;

	for(i=0; i<nCols; ++i)
	{
		PROTECT(wQvec = allocVector(INTSXP, nwRows));
		++nProtected;
		SET_VECTOR_ELT(dflist, i, wQvec);

		rQvec = VECTOR_ELT(pDf, i);
		if(TYPEOF(rQvec) != REALSXP)
			error("[get_col_quantiles] All columns in pDf must be REAL!");

		val = 0;
		k = 0; // which quantile we look for
		for(j=0; (j<nrRows) & (k<nwRows); ++j)
		{
			val += REAL(rQvec)[j];
			if(val > q[k])
			{
				// The k-th entry gets the row index
				// where the cumulative value exeeds
				// the k-th quantile
				INTEGER(wQvec)[k] = (int) j;
				++k;
			}
		}
		// Fill the rest of quantiles with max values
		for(; k<nwRows; ++k)
			INTEGER(wQvec)[k] = (int)nrRows;
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Column names
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	PROTECT(col_names = allocVector(STRSXP, nCols));
	++nProtected;
	in_col_names = getAttrib(pDf, R_NamesSymbol);

	buf = (char*) R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nCols; ++i)
	{
    	sprintf(buf, "%i", i);
    	SET_STRING_ELT(col_names, i, mkChar(CHAR(STRING_ELT(in_col_names, i))));
	}
	setAttrib(dflist, R_NamesSymbol, col_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Row names
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
    PROTECT(row_names = allocVector(STRSXP, nwRows));
    ++nProtected;

    for(i=0; i<nwRows; ++i)
    {
    	sprintf(buf, "q_%i", (int)(q[i]* 100));
    	SET_STRING_ELT(row_names, i, mkChar(buf));
    }
    setAttrib(dflist, R_RowNamesSymbol, row_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	setAttrib(dflist, R_ClassSymbol, mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP melt_vector(SEXP pValues, SEXP pFactor)
{
	int *val;
	unsigned i, j, k, n, fac, res_len;
	int sum;
	SEXP res;


	if(TYPEOF(pValues) != INTSXP)
		error("[melt_vector] pValues must be INT!");
	if(TYPEOF(pFactor) != INTSXP)
		error("[melt_vector] pFactor must be INT!");

	val = INTEGER(pValues);

	n = (unsigned) LENGTH(pValues);
	if(INTEGER(pFactor)[0] < 2)
		error("[melt_vector] Factor must be >1!");

	if( n%((unsigned)(INTEGER(pFactor)[0])) != 0)
		error("[melt_vector] Length of pValues must be multiple of %i\n", INTEGER(pFactor)[0]);

	fac = (unsigned)INTEGER(pFactor)[0];


	res_len = (unsigned) (n / fac);
	res = PROTECT(allocVector(INTSXP, res_len));
	for(i=0, k=0; i<res_len; ++i)
	{
		sum = 0;
		for(j=0; j<fac; ++j, ++k)
			sum += val[k];
		INTEGER(res)[i] = sum;
	}
	UNPROTECT(1);
	return res;
}

SEXP melt_kmer_matrix(SEXP pKmerCount, SEXP pK)
{
	int *kCount, *array;
	int old_k, new_k;
	int i, j, p, val, m;
	int nKmer_rows, nRows, nCols;
	char *buf;
	SEXP pCount_dim, pArray, dim, col_names, mat_dim_names;

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Intended for use with kmerCount matrix
	// Integer valued.
	// Rows		= DNA k mers
	// Columns	= fastq files
	//
	// Melting means that matrix is reduced in size
	// so that the result is a kmerCount matrix for lower k:
	// E.g. Counts for AAA AAA AAA ,..., AAA AAA TTT are
	// summed up into row AAA AAA.
	//
	// This shoule (presumably) reduce calculation times
	// for pairwise distances.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	if(TYPEOF(pKmerCount) != INTSXP)
		error("[melt_kmer_matrix] Matrix must be Integer!");

	if(TYPEOF(pK) != INTSXP)
		error("[melt_kmer_matrix] pK must be Integer!");

	if(length(pK) != 2)
		error("[melt_kmer_matrix] pK must contain two values!");

	kCount = INTEGER(pKmerCount);
	old_k = INTEGER(pK)[0];
	new_k = INTEGER(pK)[1];

	if(old_k <= new_k)
		error("[melt_kmer_matrix] old_k (%i) must be > new_k (%i)!", old_k, new_k);

	if((old_k > max_k) || (new_k < 1))
		error("[melt_kmer_matrix] k values must be in 1, ..., %u!", max_k);

	pCount_dim = getAttrib(pKmerCount, R_DimSymbol);
	nKmer_rows = INTEGER(pCount_dim)[0];

	if(nKmer_rows != (1 << (old_k << 1 )) )
		error("[melt_kmer_matrix] Matrix must have 4^k rows!");

	// + + + + + + + + + + + + + + + + + + + //
	// Create matrix object for returning
	// + + + + + + + + + + + + + + + + + + + //
	nRows = (1 << (new_k << 1));			// 4^new_k
	nCols = INTEGER(pCount_dim)[1];	// nFiles

	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));
	array = INTEGER(pArray);

	// + + + + + + + + + + + + + + + + + + + //
	// Dim attribute
	// + + + + + + + + + + + + + + + + + + + //
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	// + + + + + + + + + + + + + + + + + + + //
	// Column names: 1:nCols
	// + + + + + + + + + + + + + + + + + + + //
	col_names = PROTECT(allocVector(STRSXP, nCols));
	buf = R_alloc(rbuf_size, sizeof(char));
	for(i=0; i<nCols; ++i)
	{
		sprintf(buf, "%i", i + 1);
		SET_STRING_ELT(col_names, i, mkChar(buf));
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Row names: DNA k-mers
	// + + + + + + + + + + + + + + + + + + + //
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(mat_dim_names, 0, PROTECT(create_dna_k_mers(new_k)));
	SET_VECTOR_ELT(mat_dim_names, 1, col_names);
	setAttrib(pArray, R_DimNamesSymbol, mat_dim_names);

	// + + + + + + + + + + + + + + + + + + + //
	// Do melt down
	// + + + + + + + + + + + + + + + + + + + //
	m = (1 << ((old_k - new_k) << 1));
	for(i=0, p=0; i < nRows * nCols; ++i)
	{
		val = 0;
		for(j=0; j<m; ++j, ++p)
			val += kCount[p];
		array[i] = val;
	}
	UNPROTECT(5);
	return pArray;
}

SEXP scale_kmer_matrix(SEXP pKmerCount, SEXP pScale)
{
	int nRows, nCols;
	int i, column_index, row_index;
	int *out, *in;
	double scale;

	SEXP pInDim, dim, pArray;
	SEXP pInDimNames, pInRowNames, pInColNames;
	SEXP pOutDimNames, pOutRowNames, pOutColNames;

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Intended for use with kmerCount matrix
	// Integer valued.
	// Rows		= DNA k mers
	// Columns	= fastq files
	//
	// Scaling the matrix intends that count values
	// in each column are multiplicated
	// (by a value >1 i.e. scaled up)
	// so that the column sums become nearly equal.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	if(TYPEOF(pKmerCount) != INTSXP)
		error("[scale_kmer_matrix] Matrix must be Integer!");

	pInDim = getAttrib(pKmerCount, R_DimSymbol);
	nRows = INTEGER(pInDim)[0];
	nCols = INTEGER(pInDim)[1];

	if(TYPEOF(pScale) != REALSXP)
		error("[scale_kmer_matrix] pScale must be Real!");

	if(length(pScale) != nCols)
		error("[scale_kmer_matrix] There must be one scale for each matrix column!");

	// + + + + + + + + + + + + + + + + + + + //
	// Create output matrix
	// + + + + + + + + + + + + + + + + + + + //
	pArray = PROTECT(allocMatrix(INTSXP, nRows, nCols));

	// + + + + + + + + + + + + + + + + + + + //
	// Insert values into output matrix
	// + + + + + + + + + + + + + + + + + + + //
	column_index = 0;
	row_index = 0;
	out = INTEGER(pArray);
	in = INTEGER(pKmerCount);

	// + + + + + + + + + + + + + + + + + + + //
	// Initialize scale
	// + + + + + + + + + + + + + + + + + + + //
	scale = REAL(pScale)[column_index];
	if(scale < 1)
		error("[scale_kmer_matrix] scale[%i]=%i must be >1!", column_index, scale);

	for(i=0; i < (nRows * nCols); ++i)
	{
		if(row_index == nRows)
		{
			row_index = 0;
			++column_index;
			scale = REAL(pScale)[column_index];
			if(scale < 1)
				error("[scale_kmer_matrix] scale[%i] = %i must be >1!", column_index, scale);
		}
		out[i] = (int)(in[i] * scale);
		++row_index;
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Dim symbol
	// + + + + + + + + + + + + + + + + + + + //
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = nRows;
	INTEGER(dim)[1] = nCols;
	setAttrib(pArray, R_DimSymbol, dim);

	// + + + + + + + + + + + + + + + + + + + //
	// Copy dim names
	// + + + + + + + + + + + + + + + + + + + //
	pInDimNames = getAttrib(pKmerCount, R_DimNamesSymbol);

	pInRowNames = VECTOR_ELT(pInDimNames, 0);
	pOutRowNames = PROTECT(allocVector(STRSXP, nRows));

	for(i=0;i<nRows;++i)
		SET_STRING_ELT(pOutRowNames, i, mkChar(CHAR(STRING_ELT(pInRowNames, i))));

	pInColNames = VECTOR_ELT(pInDimNames, 1);
	pOutColNames = PROTECT(allocVector(STRSXP, nCols));
	for(i=0; i<nCols; ++i)
		SET_STRING_ELT(pOutColNames, i, mkChar(CHAR(STRING_ELT(pInColNames, i))));

	pOutDimNames = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(pOutDimNames, 0, pOutRowNames);
	SET_VECTOR_ELT(pOutDimNames, 1, pOutColNames);
	setAttrib(pArray, R_DimNamesSymbol, pOutDimNames);


	UNPROTECT(5);
	return pArray;
}


static R_INLINE double cb_distance(int * x, int * y, int size)
{
	int i;
	double nom, denom;
	double val;
	val = 0;
	for(i=0; i<size; ++i)
	{
		if( (x[i] < 0) || (y[i] < 0) )
			error("[cb_distance] No negative values allowed!");
		if(!( (x[i] > 0) || (y[i] > 0)))
			return 0;

		nom	 = (double)( (x[i] > y[i]) ? (x[i]-y[i]) : (y[i]-x[i]) );
		denom= (double)(x[i] + y[i]);
		val += nom / denom;
	}
	return val / size;
}

SEXP cb_distance_matrix(SEXP pKmerCount)
{
	// ToDo: Check for matrix??
	int i, j, size, nKmer_rows;
	int *kCount;
	double *array;
	size_t array_size;

	SEXP count_dim, res, dim, mat_dim_names;
	SEXP pInDimNames, pInColNames, dim_names;


	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Intended for use with kmerCount matrix
	// Integer valued.
	// Rows		= DNA k mers
	// Columns	= fastq files
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //

	if(TYPEOF(pKmerCount) != INTSXP)
		error("[cb_distance_matrix] Count matrix must be Int!");

	count_dim = getAttrib(pKmerCount, R_DimSymbol);
	kCount = INTEGER(pKmerCount);

	size = INTEGER(count_dim)[1];
	nKmer_rows = INTEGER(count_dim)[0];

	// + + + + + + + + + + + + + + + + + + + //
	// Create matrix for returning
	// + + + + + + + + + + + + + + + + + + + //
	res = PROTECT(allocMatrix(REALSXP, size, size));
	array = REAL(res);
	array_size = ((size_t)(size * size));
	memset(array, 0, array_size * sizeof(double));

	// + + + + + + + + + + + + + + + + + + + //
	// Dim attribute
	// + + + + + + + + + + + + + + + + + + + //
	dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = size;
	INTEGER(dim)[1] = size;
	setAttrib(res, R_DimSymbol, dim);

	// + + + + + + + + + + + + + + + + + + + //
	// Dim names
	// + + + + + + + + + + + + + + + + + + + //
	mat_dim_names = PROTECT(allocVector(VECSXP, 2));

	// + + + + + + + + + + + + + + + + + + + //
	// Copy dim names
	// + + + + + + + + + + + + + + + + + + + //
	pInDimNames = getAttrib(pKmerCount, R_DimNamesSymbol);
	pInColNames = VECTOR_ELT(pInDimNames, 1);

	// + + + + + + + + + + + + + + + + + + + //
	// Row names
	// + + + + + + + + + + + + + + + + + + + //
	PROTECT(dim_names = allocVector(STRSXP, size));
	for(i=0; i<size; ++i)
		SET_STRING_ELT(dim_names, i, mkChar(CHAR(STRING_ELT(pInColNames, i))));
	SET_VECTOR_ELT(mat_dim_names, 0, dim_names);

	// + + + + + + + + + + + + + + + + + + + //
	// Col names
	// + + + + + + + + + + + + + + + + + + + //
	PROTECT(dim_names = allocVector(STRSXP, size));
	for(i=0; i<size; ++i)
		SET_STRING_ELT(dim_names, i, mkChar(CHAR(STRING_ELT(pInColNames, i))));

	SET_VECTOR_ELT(mat_dim_names, 1, dim_names);
	setAttrib(res, R_DimNamesSymbol, mat_dim_names);

	// + + + + + + + + + + + + + + + + + + + //
	// Iterate through res matrix and
	// insert distances
	// + + + + + + + + + + + + + + + + + + + //

	for(i=0; i<size; ++i)			// i = column index
	{
		for(j=i+1; j<size; ++j)	// j= row index
			array[(i * size) + j] = cb_distance(kCount + (nKmer_rows * j), kCount + (nKmer_rows * i), nKmer_rows);
	}
	UNPROTECT(5);
	return res;
}

SEXP rel_int(SEXP pInt)
{
	unsigned i, n;
	unsigned long int sum;
	int *val;

	SEXP res;

	if(TYPEOF(pInt) != INTSXP)
		error("[rel_int] pInt must be Integer!");

	n = (unsigned)length(pInt);
	val = INTEGER(pInt);
	sum = 0;
	for(i=0; i<n; ++i)
		sum += (unsigned long) val[i];

	if(!sum)
		return pInt;

	res = PROTECT(allocVector(REALSXP, n));
	for(i=0; i<n; ++i)
		REAL(res)[i] = ((double)val[i]/((double)sum));

	UNPROTECT(1);
	return res;
}

SEXP rel_real(SEXP pReal)
{
	unsigned i, n;
	double sum;
	double *val;
	SEXP res;

	if(TYPEOF(pReal) != REALSXP)
		error("[rel_real] pReal must be REAL!");

	n = (unsigned)length(pReal);
	val = REAL(pReal);
	sum = 0;
	for(i=0; i<n; ++i)
		sum += val[i];
	if(!sum)
		return pReal;

	res = PROTECT(allocVector(REALSXP, n));
	for(i=0; i<n; ++i)
		REAL(res)[i] = (val[i] / sum);

	UNPROTECT(1);
	return res;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 	Simulation of fastq Reads
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

SEXP sim_k_values(SEXP pVal)
{
	int k, n;
	double max_val;
	SEXP res;

	if(TYPEOF(pVal) != INTSXP)
		error("[sim_k_values] pVal must be Int!");

	if(LENGTH(pVal) != 2)
		error(" [sim_k_values] pVal must have length 2!");

	k = INTEGER(pVal)[0];
	n = INTEGER(pVal)[1];
	max_val = (double) (1 << (k << 1));

	res = PROTECT(allocVector(REALSXP, n));
	rsim_unif(n, 0, max_val, REAL(res));

	UNPROTECT(1);
	return res;
}

SEXP sim_dna_k_mer(SEXP pVal)
{
	int k, nk, nSeq;
	char *buf, *buf_iter;
	double max_kmer_index;
	int kmer_index;
	int i, j, p, q;
	size_t seqlen;
	SEXP res;

	if(TYPEOF(pVal) != INTSXP)
		error("[sim_k_values] pVal must be Int!");

	if(LENGTH(pVal) != 3)
		error(" [sim_k_values] pVal must have length 3!");

	k = INTEGER(pVal)[0];
	nk = INTEGER(pVal)[1];
	nSeq = INTEGER(pVal)[2];
	seqlen = (size_t) (k * nk);

	res = PROTECT(allocVector(STRSXP, nSeq));

	buf = R_alloc(seqlen+1, sizeof(char));
	buf[seqlen] = '\0';
	max_kmer_index = (double) (1 << (k << 1));

	GetRNGstate();
	for(i=0; i<nSeq; ++i)
	{
		for(j=0; j<nk; ++j)
		{
			// put j'th kmer into buffer
			buf_iter = buf + j * k;
			kmer_index = (int)runif(0, max_kmer_index);
			//Rprintf("[sim_dna_k_mer] kmer_index=%i\n", kmer_index);

			q = (2 * (k - 1));
			for(p=0; p<k; ++p)
			{
				buf_iter[p] = (char) rev_ACGT[(kmer_index >> q) & 3];
				q -= 2;
			}
			//Rprintf("[sim_dna_k_mer] buf_iter='%s'\n", buf_iter);
		}
		SET_STRING_ELT(res, i, mkChar(buf));
	}
	PutRNGstate();
	UNPROTECT(1);
	return res;
}

SEXP set_dna_k_mer(SEXP pVal, SEXP pSeq, SEXP pPosition, SEXP pKmerIndex, SEXP pSeqIndex)
{
	int k, nk, nSeq, nPos, seqlen, nInsSeq;
	int *position, *kmer_index;
	int i, j, p, q;
	char *buf, *buf_iter;
	const char *seq;
	SEXP res;

	if(TYPEOF(pVal) != INTSXP)
		error("[set_dna_k_mer] pVal must be Int!");

	if(LENGTH(pVal) != 3)
		error(" [set_dna_k_mer] pVal must have length 3!");

	if(TYPEOF(pSeq) != STRSXP)
		error("[set_dna_k_mer] pSeq must be String!");

	if(TYPEOF(pPosition) != INTSXP)
		error("[set_dna_k_mer] pPosition must be Int!");

	if(TYPEOF(pKmerIndex) != INTSXP)
		error("[set_dna_k_mer] pKmerIndex must be Int!");

	if(TYPEOF(pSeqIndex) !=  INTSXP)
		error("[set_dna_k_mer] pSeqIndex must be Int!");


	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Sequences consist of nk k-mer motifs.
	// They make up the random part of motifs
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	k = INTEGER(pVal)[0];
	nk = INTEGER(pVal)[1];
	seqlen = k * nk;
	nSeq = LENGTH(pSeq);


	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// At 'position' a new k-mer is inserted
	// This will be the 'deterministic' part of motifs.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	nPos = LENGTH(pPosition);
	if(LENGTH(pKmerIndex) != nPos)
		error("[set_dna_k_mer] pPosition and pKmerIndex must have equal length!");

	position = INTEGER(pPosition);
	kmer_index = INTEGER(pKmerIndex);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// The deterministic k-mers will be only inserted
	// in the first sequences
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	nInsSeq = INTEGER(pSeqIndex)[0];
	if(nInsSeq<1 || nInsSeq>nSeq)
		error("[seq_dna_k_mer] pSeqIndex out of range!");

	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// The function returns a modified copy of the
	// given sequence data.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + //
	res = PROTECT(allocVector(STRSXP, nSeq));
	buf = R_alloc(((size_t)seqlen)+1, sizeof(char));
	buf[seqlen] = '\0';

	for(i=0; i<nSeq; ++i)
	{
		for(j=0; j<nPos; ++j)
		{
			if( (position[j] < 0) || (position[j] > seqlen - k + 1) )
				error("[set_dna_k_mer] position[%i] out of range!", j);

			// copy sequence
			seq = CHAR(STRING_ELT(pSeq, i));
			for(p=0; p<seqlen; ++p)
				buf[p] = seq[p];

			if(i < nInsSeq)
			{
				buf_iter = buf + position[j];
				q = (2 * (k - 1));
				for(p=0; p<k; ++p)
				{
					buf_iter[p] = (char) rev_ACGT[ (kmer_index[j] >> q) & 3];
					q -= 2;
				}
			}
		}
		SET_STRING_ELT(res, i, mkChar(buf));
	}
	UNPROTECT(1);
	return res;
}


SEXP gzwrite_fastq_dna(SEXP pVal, SEXP pSeq, SEXP pFilename)
{
	int i, k, n, nk, nBytes;
	char *qSeq, *headBuf;
	size_t seqlen;

	gzFile gz;
	SEXP res;

	if(TYPEOF(pVal) != INTSXP)
		error("[gzwrite_fastq_dna] pVal must be Int!");

	if(LENGTH(pVal) != 3)
		error("[gzwrite_fastq_dna] pVal must have length 3!");

	if(TYPEOF(pFilename) != STRSXP)
		error("[gzwrite_fastq_dna] pFilename must be string!");

	if(TYPEOF(pSeq) != STRSXP)
		error("[gzwrite_fastq_dna] pSeq must be string!");

	k = INTEGER(pVal)[0];		    // pval[0] = k
	nk = INTEGER(pVal)[1];	        // pval[1] = nk
	//int nSeq=INTEGER(pVal)[2];	// pval[2] = nSeq
	seqlen = (size_t) (k * nk);

	// All phred values will be '.'=46-33=13
	qSeq = R_alloc(seqlen + 1, sizeof(char));
	memset(qSeq, 46, seqlen* sizeof(char));
	qSeq[seqlen] = '\0';

	gz = gzopen(CHAR(STRING_ELT(pFilename, 0)), "wb");
	if(!gz)
		error("[gzwrite_fastq_dna] Could not open file '%s'!", CHAR(STRING_ELT(pFilename, 0)));

	nBytes = 0;
	n = LENGTH(pSeq);
	headBuf = R_alloc(50, sizeof(char));

	for(i=0; i<n; ++i)
	{
		sprintf(headBuf, "@%i\n", i);
		nBytes += gzputs(gz, headBuf);
		nBytes += gzputs(gz, CHAR(STRING_ELT(pSeq, i)));
		nBytes += gzputs(gz, "\n+\n");
		nBytes += gzputs(gz, qSeq);
		nBytes += gzputs(gz, "\n");
	}
	gzclose(gz);

	res=PROTECT(allocVector(INTSXP, 1));
	INTEGER(res)[0] = nBytes;
	UNPROTECT(1);
	return res;
}

SEXP gzwrite_mod_fastq_dna(SEXP pVal, SEXP pFilename, SEXP pSeq)
{
	int i, n, k, nk, nBytes;
	char *qSeq, *headBuf;
	gzFile gz;
	SEXP res;
	size_t seqlen;

	if(TYPEOF(pVal) != INTSXP)
		error("[gzwrite_fastq_dna] pVal must be Int!");

	if(LENGTH(pVal) != 3)
		error("[gzwrite_fastq_dna] pVal must have length 3!");

	if(TYPEOF(pFilename) != STRSXP)
		error("[gzwrite_fastq_dna] pFilename must be string!");

	if(TYPEOF(pSeq) != STRSXP)
		error("[gzwrite_fastq_dna] pSeq must be string!");

	k = INTEGER(pVal)[0];		 // pval[0] = k
	nk=INTEGER(pVal)[1];	     // pval[1] = nk
	//int nSeq=INTEGER(pVal)[2]; // pval[2] = nSeq
	seqlen= (size_t) (k * nk);

	// All phred values will be '.'=46-33=13
	qSeq = R_alloc(seqlen + 1, sizeof(char));
	memset(qSeq, 46, seqlen* sizeof(char));
	qSeq[seqlen] = '\0';

	gz = gzopen(CHAR(STRING_ELT(pFilename, 0)), "wb");
	if(!gz)
		error("[gzwrite_fastq_dna] Could not open file '%s'!", CHAR(STRING_ELT(pFilename, 0)));

	nBytes = 0;
	n = LENGTH(pSeq);

	headBuf = R_alloc(50, sizeof(char));
	for(i=0; i<n; ++i)
	{
		sprintf(headBuf, "@%i\n", i);
		nBytes += gzputs(gz, headBuf);
		nBytes += gzputs(gz, CHAR(STRING_ELT(pSeq, 0)));
		nBytes += gzputs(gz, "\n+\n");
		nBytes += gzputs(gz, qSeq);
		nBytes += gzputs(gz, "\n");
	}
	gzclose(gz);

	res = PROTECT(allocVector(INTSXP, 1));
	INTEGER(res)[0] = nBytes;
	UNPROTECT(1);
	return res;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 	R_init
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


void R_init_seqTools(DllInfo * info)
{
	R_CallMethodDef cmd[] ={
			{ "r_enlarge_int_mat", 		(DL_FUNC) &r_enlarge_int_mat,		2},

			{ "create_dna_k_mers",		(DL_FUNC) &create_dna_k_mers,       1},
			{ "count_Kmers",			(DL_FUNC) &count_Kmers,				4},
			{ "count_dna_Kmers",		(DL_FUNC) &count_dna_Kmers,			5},
			{ "rev_count_dna_Kmers",	(DL_FUNC) &rev_count_dna_Kmers,		5},
			{ "count_genome_Kmers",		(DL_FUNC) &count_genome_Kmers,		7},
			{ "count_splice_Kmers",		(DL_FUNC) &count_splice_Kmers,		8},
			{ "count_fastq",			(DL_FUNC) &count_fastq,				2},
			{ "trim_fastq",				(DL_FUNC) &trim_fastq,				3},
			{ "fastq_KmerSubset_locs",	(DL_FUNC) &fastq_KmerSubset_locs,	3},
			{ "fastq_Kmer_locs",		(DL_FUNC) &fastq_Kmer_locs,			2},
			{ "write_fai",				(DL_FUNC) &write_fai,				2},
			{ "count_fasta_Kmers",		(DL_FUNC) &count_fasta_Kmers,		2},

			{ "get_column_quantiles",	(DL_FUNC) &get_column_quantiles,	2},
			{ "melt_vector",			(DL_FUNC) &melt_vector,				2},
			{ "melt_kmer_matrix",		(DL_FUNC) &melt_kmer_matrix,		2},
			{ "scale_kmer_matrix",		(DL_FUNC) &scale_kmer_matrix,		2},
			{ "cb_distance_matrix",		(DL_FUNC) &cb_distance_matrix,		1},
			{ "rel_int",				(DL_FUNC) &rel_int,					1},
			{ "rel_real",				(DL_FUNC) &rel_real,				1},

			{ "sim_k_values",			(DL_FUNC) &sim_k_values,			1},
			{ "sim_dna_k_mer",			(DL_FUNC) &sim_dna_k_mer,			1},
			{ "set_dna_k_mer",			(DL_FUNC) &set_dna_k_mer,			5},
			{ "gzwrite_fastq_dna",		(DL_FUNC) &gzwrite_fastq_dna,		3},
			{NULL, NULL, 0}
	};
	//			{ "",	(DL_FUNC) &,	}
	R_registerRoutines(info, NULL, cmd, NULL, NULL);
}



#endif /* SEQTOOLS_C_ */
