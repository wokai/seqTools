/*
 * stat_defs.h
 *
 *  Created on: 14.10.2013
 *      Author: kaisers
 */

#ifndef STAT_DEFS_H_
#define STAT_DEFS_H_

/*
 * max_k is set to 15: 4^15 = 1.073.741.824 combinations = 8 GByte.
 * This allows usage of (signed) int as array index.
 */
static const int max_k = 15;

/*
 * 93 is the maximal theoretically allowed phred value.
 */
static const int maxPhred = 93;

/*
 * Buffer size for printing ints into chars
 */

// Must have same size:
static const int buf_size = 1024;			// Version for gzgets
static const size_t rbuf_size = 1024;		// Version for R_alloc

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Static variable declararions for maxEntScores (scoreconsensus and maxent_score functions)
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// Translation tables: DNA-Nuc -> integer index
# define aval 0
# define cval 1
# define gval 2
# define tval 3
# define nval 4
# define zval 5
# define lfcv 6		/* line feed */

// ASCII:	A 65, C 67, G 71, T 84, a 97, c 99, g 103, t 116
// Translation:	A->0, C->1, G->2, T->3
static const unsigned char ACGT[256] = {
	// 64 - block																					// start index
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 0
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 16
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 32
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 48

	zval,aval,zval,cval,	zval,zval,zval,gval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 64
	zval,zval,zval,zval,	tval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 80
	zval,aval,zval,cval,	zval,zval,zval,gval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 96
	zval,zval,zval,zval,	tval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 112

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 128
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 144
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 160
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 176

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 192
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 208
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 224
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval		// 240
};

// ASCII:	A 65, C 67, G 71, T 84, a 97, c 99, g 103, t 116
// Translation:	A->0, C->1, G->2, T->3
static const unsigned char ACGTN[256] = {
	// 64 - block																					// start index
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 0
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 16
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 32
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 48

	zval,aval,zval,cval,	zval,zval,zval,gval,	zval,zval,zval,zval,	zval,zval,nval,zval,	// 64
	zval,zval,zval,zval,	tval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 80
	zval,aval,zval,cval,	zval,zval,zval,gval,	zval,zval,zval,zval,	zval,zval,nval,zval,	// 96
	zval,zval,zval,zval,	tval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 112

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 128
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 144
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 160
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 176

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 192
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 208
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 224
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval		// 240
};


# define zzv 0
# define gcv 1

static const unsigned char GCZ[256] = {
	// 64 - block																	// start index
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 0
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 16
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 32
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 48

	zzv,zzv,zzv,gcv,	zzv,zzv,zzv,gcv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 64
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 80
	zzv,zzv,zzv,gcv,	zzv,zzv,zzv,gcv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 96
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 112

	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 128
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 144
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 160
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 176

	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 192
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 208
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	// 224
	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv,	zzv,zzv,zzv,zzv		// 240
};




// Translation into reverse sequence values
// ASCII:	A 65, C 67, G 71, T 84, a 97, c 99, g 103, t 116
// Translation:	A->T->3, C->G->2, G->C->1, T->A->0
static const unsigned char revACGT[256] = {
	// 64 - block																		// start index
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 0
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 16
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 32
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 48

	zval,tval,zval,gval,	zval,zval,zval,cval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 64
	zval,zval,zval,zval,	aval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 80
	zval,tval,zval,gval,	zval,zval,zval,cval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 96
	zval,zval,zval,zval,	aval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 112

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 128
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 144
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 160
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 176

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 192
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 208
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 224
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval		// 240
};


// Character encoding for create_dna_k_mers
static const unsigned char rev_ACGT[4] = { 65, 67, 71, 84 };

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// IUPAC encodings used for fastq
//
///////////////////////////////////////////////////////////////////////////////////////////////////

//	Letter		Base				Code
//	A			Adenine				0
//	C			Cytosine			1
//	G			Guanine				2
//	T			Thymine				3
//	U			Uracil				4
//	R			A or G				5
//	Y			C or T				6
//	S			G or C				7
//	W			A or T				8
//	K			G or T				9
//	M			A or C				10
//	B			C or G or T			11
//	D			A or G or T			12
//	H			A or C or T			13
//	V			A or C or G			14
//	N			any base			15
//	.			Gap					16
//	-			Gap					17

//	=			Gap					18
//  Any other value					20


# define iua 0
# define iuc 1
# define iug 2
# define iut 3
# define iuu 4
# define iur 5
# define iuy 6
# define ius 7
# define iuw 8
# define iuk 9
# define ium 10
# define iub 11
# define iud 12
# define iuh 13
# define iuv 14
# define iun 15
# define iup 16
# define iux 17
# define iue 18
# define iuz 20

static int const nIupac = 19;
static const char * iupac_chars = "acgturyswkmbdhvn.-=";


static const unsigned char IUPAC[256] = {
	// 64 - block																	// start index
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 0
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 16
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iux,iup,iuz,	// 32
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iue,iuz,iuz,	// 48

	iuz,iua,iub,iuc,	iud,iuz,iuz,iug,	iuh,iuz,iuz,iuk,	iuz,ium,iun,iuz,	// 64
	iuz,iuz,iur,ius,	iut,iuu,iuv,iuw,	iux,iuy,iuz,iuz,	iuz,iuz,iuz,iuz,	// 80
	iuz,iua,iub,iuc,	iud,iuz,iuz,iug,	iuh,iuz,iuz,iuk,	iuz,ium,iun,iuz,	// 96
	iuz,iuz,iur,ius,	iut,iuu,iuv,iuw,	iux,iuy,iuz,iuz,	iuz,iuz,iuz,iuz,	// 112

	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 128
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 144
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 160
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 176

	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 192
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 208
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	// 224
	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz,	iuz,iuz,iuz,iuz		// 240
};

static const unsigned char LFCR[256] = {
	// 64 - block																					// start index
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,lfcv,zval,	zval,lfcv,zval,zval,	// 0
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 16
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 32
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 48

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 64
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 80
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 96
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 112

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 128
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 144
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 160
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 176

	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 192
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 208
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	// 224
	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval,	zval,zval,zval,zval		// 240
};



#endif /* STAT_DEFS_H_ */
