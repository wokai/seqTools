/*
 * resize_matrix.h
 *
 *  Created on: 05.11.2013
 *      Author: wolfgang
 */

#ifndef RESIZE_MATRIX_H_
#define RESIZE_MATRIX_H_

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>

SEXP r_enlarge_int_mat(SEXP pMat, SEXP pNewDim)
{
	if(TYPEOF(pMat) != INTSXP)
		error("[r_enlarge_int_mat] pMat matrix must be integer!");

	if(TYPEOF(pNewDim) != INTSXP)
		error("[r_enlarge_int_mat] pNewDim must be integer!");

	if(LENGTH(pNewDim) != 2)
		error("[r_enlarge_int_mat] pNewDim must have length 2!");

	int new_nRows = INTEGER(pNewDim)[0];
	int new_nCols = INTEGER(pNewDim)[1];

	if( (new_nRows < 1) || (new_nCols < 1) )
		error("[r_enlarge_int_mat] pNewDim values must be positive!");

	size_t array_size = (size_t) (new_nRows * new_nCols);
	SEXP pMat_dim = getAttrib(pMat, R_DimSymbol);
	int *mat = INTEGER(pMat);

	int old_nRows = INTEGER(pMat_dim)[0];
	int old_nCols = INTEGER(pMat_dim)[1];

	// + + + + + + + + + + + + + + + + + + + //
	// Do nothing when new dim's are not larger
	// + + + + + + + + + + + + + + + + + + + //
	if( (new_nRows < old_nRows) || (new_nCols < old_nCols) )
		return pMat;

	if( (new_nRows == old_nRows) && (new_nCols == old_nCols) )
		return pMat;

	SEXP pRes = PROTECT(allocVector(INTSXP, (R_xlen_t) array_size));

	int *res = INTEGER(pRes);
	memset(res, 0, array_size * sizeof(int));


	int i, j, old_offset, new_offset;
	for(i=0; i<old_nCols; ++i)
	{
		old_offset = i * old_nRows;
		new_offset = i * new_nRows;

		for(j=0; j<old_nRows; ++j)
			res[ new_offset + j ] = mat[ old_offset + j ];
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Dim attribute
	// + + + + + + + + + + + + + + + + + + + //
	SEXP dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = new_nRows;
	INTEGER(dim)[1] = new_nCols;
	setAttrib(pRes, R_DimSymbol, dim);


	// + + + + + + + + + + + + + + + + + + + //
	// Copy dim names
	// + + + + + + + + + + + + + + + + + + + //
	char *buf = R_alloc(rbuf_size, sizeof(char));
	SEXP pInDimNames = getAttrib(pMat, R_DimNamesSymbol);

	SEXP pInRowNames = VECTOR_ELT(pInDimNames, 0);
	SEXP pOutRowNames = PROTECT(allocVector(STRSXP, new_nRows));
	for(i=0; i<old_nRows; ++i)
		SET_STRING_ELT(pOutRowNames, i, mkChar(CHAR(STRING_ELT(pInRowNames, i))));
	for(i=old_nRows; i<new_nRows; ++i)
	{
		sprintf(buf, "%i", i + 1);
		SET_STRING_ELT(pOutRowNames, i, mkChar(buf));
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Column names: 1:nCols
	// + + + + + + + + + + + + + + + + + + + //
	SEXP pInColNames=VECTOR_ELT(pInDimNames, 1);
	SEXP pOutColNames=PROTECT(allocVector(STRSXP, new_nCols));

	for(i=0; i<old_nCols; ++i)
		SET_STRING_ELT(pOutColNames, i, mkChar(CHAR(STRING_ELT(pInColNames, i))));

	for(i=old_nCols; i<new_nCols; ++i)
	{
		sprintf(buf, "%i", i+1);
		SET_STRING_ELT(pOutColNames, i, mkChar(buf));
	}

	SEXP pOutDimNames=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(pOutDimNames, 0, pOutRowNames);
	SET_VECTOR_ELT(pOutDimNames, 1, pOutColNames);
	setAttrib(pRes, R_DimNamesSymbol, pOutDimNames);

	UNPROTECT(5);
	return pRes;
}



SEXP enlarge_int_mat(SEXP pMat, int new_nRows, int new_nCols)
{
	if(TYPEOF(pMat)!=INTSXP)
		error("[enlarge_int_mat] pMat matrix must be integer!");

	size_t array_size = (size_t) (new_nRows * new_nCols);
	SEXP pMat_dim = getAttrib(pMat, R_DimSymbol);

	int *mat = INTEGER(pMat);
	int old_nRows = INTEGER(pMat_dim)[0];
	int old_nCols = INTEGER(pMat_dim)[1];

	if( (new_nRows < old_nRows) || (new_nCols < old_nCols) )
		return pMat;

	if( (new_nRows == old_nRows) && (new_nCols == old_nCols) )
		return pMat;

	SEXP pRes = PROTECT(allocVector(INTSXP, (R_xlen_t)array_size));

	int *res = INTEGER(pRes);
	memset(res, 0, array_size * sizeof(int));

	int i, j, old_offset, new_offset;
	for(i=0; i<old_nCols; ++i)
	{
		old_offset = i * old_nRows;
		new_offset = i * new_nRows;

		for(j=0; j<old_nRows; ++j)
			res[ new_offset + j ] = mat[ old_offset + j ];
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Dim attribute
	// + + + + + + + + + + + + + + + + + + + //
	SEXP dim = PROTECT(allocVector(INTSXP, 2));

	INTEGER(dim)[0] = new_nRows;
	INTEGER(dim)[1] = new_nCols;

	setAttrib(pRes, R_DimSymbol, dim);


	// + + + + + + + + + + + + + + + + + + + //
	// Copy dim names
	// + + + + + + + + + + + + + + + + + + + //
	char *buf = R_alloc(rbuf_size, sizeof(char));
	SEXP pInDimNames = getAttrib(pMat, R_DimNamesSymbol);

	SEXP pInRowNames = VECTOR_ELT(pInDimNames, 0);
	SEXP pOutRowNames = PROTECT(allocVector(STRSXP, new_nRows));

	for(i=0; i<old_nRows; ++i)
		SET_STRING_ELT(pOutRowNames, i, mkChar(CHAR(STRING_ELT(pInRowNames, i))));

	for(i=old_nRows; i<new_nRows; ++i)
	{
		sprintf(buf, "%i", i+1);
		SET_STRING_ELT(pOutRowNames, i, mkChar(buf));
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Column names: 1:nCols
	// + + + + + + + + + + + + + + + + + + + //
	SEXP pInColNames = VECTOR_ELT(pInDimNames, 1);

	SEXP pOutColNames = PROTECT(allocVector(STRSXP, new_nCols));
	for(i=0; i<old_nCols; ++i)
		SET_STRING_ELT(pOutColNames, i, mkChar(CHAR(STRING_ELT(pInColNames, i))));

	for(i=old_nCols; i<new_nCols; ++i)
	{
		sprintf(buf, "%i", i+1);
		SET_STRING_ELT(pOutColNames, i, mkChar(buf));
	}

	SEXP pOutDimNames=PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(pOutDimNames, 0, pOutRowNames);
	SET_VECTOR_ELT(pOutDimNames, 1, pOutColNames);
	setAttrib(pRes, R_DimNamesSymbol, pOutDimNames);

	UNPROTECT(5);
	return pRes;
}

SEXP cut_down_int_mat(SEXP pMat, int new_nRows, int new_nCols)
{
	if(TYPEOF(pMat)!=INTSXP)
		error("[cut_down_int_mat] pMat matrix must be integer!");

	int array_size = new_nRows * new_nCols;

	SEXP pMat_dim = getAttrib(pMat, R_DimSymbol);

	int *mat = INTEGER(pMat);

	int old_nRows = INTEGER(pMat_dim)[0];
	int old_nCols = INTEGER(pMat_dim)[1];

	if( (new_nRows > old_nRows) || (new_nCols > old_nCols) )
		return pMat;

	if( (new_nRows == old_nRows) && (new_nCols == old_nCols) )
		return pMat;

	SEXP pRes = PROTECT(allocVector(INTSXP, array_size));
	int *res = INTEGER(pRes);


	int i, j, old_offset, new_offset;
	for(i=0; i<new_nCols; ++i)
	{
		old_offset = i * old_nRows;
		new_offset = i * new_nRows;

		for(j=0; j<new_nRows; ++j)
			res[ new_offset + j ] = mat[ old_offset + j ];
	}

	// + + + + + + + + + + + + + + + + + + + //
	// Dim attribute
	// + + + + + + + + + + + + + + + + + + + //
	SEXP dim = PROTECT(allocVector(INTSXP, 2));
	INTEGER(dim)[0] = new_nRows;
	INTEGER(dim)[1] = new_nCols;
	setAttrib(pRes, R_DimSymbol, dim);

	// + + + + + + + + + + + + + + + + + + + //
	// Copy dim names
	// + + + + + + + + + + + + + + + + + + + //
	SEXP pInDimNames = getAttrib(pMat, R_DimNamesSymbol);

	SEXP pInRowNames = VECTOR_ELT(pInDimNames, 0);
	SEXP pOutRowNames = PROTECT(allocVector(STRSXP, new_nRows));

	for(i=0; i<new_nRows; ++i)
		SET_STRING_ELT(pOutRowNames, i, mkChar(CHAR(STRING_ELT(pInRowNames, i))));

	SEXP pInColNames=VECTOR_ELT(pInDimNames, 1);
	SEXP pOutColNames=PROTECT(allocVector(STRSXP, new_nCols));

	for(i=0; i<new_nCols; ++i)
		SET_STRING_ELT(pOutColNames, i, mkChar(CHAR(STRING_ELT(pInColNames, i))));

	SEXP pOutDimNames = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(pOutDimNames, 0, pOutRowNames);
	SET_VECTOR_ELT(pOutDimNames, 1, pOutColNames);
	setAttrib(pRes, R_DimNamesSymbol, pOutDimNames);

	UNPROTECT(5);
	return pRes;
}

#endif /* RESIZE_MATRIX_H_ */
