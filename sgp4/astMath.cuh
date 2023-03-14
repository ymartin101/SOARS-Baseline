#ifndef _astMath_h_
#define _astMath_h_
/* --------------------------------------------------------------------
*
*                                astMath.cuh
*
*    this file contains miscallaneous math functions. vectors use the usual 0-1-2
*    indexing scheme.
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              11 dec 07  david vallado
*                           fixes to matrix operations
*              10 aug 06  david vallado
*                           use fmod
*              20 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
* ----------------------------------------------------------------------      */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

#pragma once

namespace astMath
{
	//global interfaces
#define pi 3.14159265358979323846
#define infinite  999999.9
#define undefined 999999.1

	/*
	*  because c++ has no multi-dimensioned, dynamic allocated matrix, I built one.
	*  the matrix type uses the standard vector to allow dynamic re-sizing during runtime.
	*  the first time a variable is defined, you define it as follows
	*       matrix(4,3);
	*  afterwards, you can just pass the variable.
	*       matrix mat3;
	*  anytime you create a new matrix inside a procedure or function, you must re-size
	*  the matrix. use the following statements to re-size!!
	*       mat3.resize(rows);  // set the # of rows
	*       for (std::vector< std::vector<double> >::iterator it=mat3.begin(); it != mat3.end();++it)
	*           it->resize(cols);  // set the # of cols
	*  the type definition has no size so you can use any size you like ... and then change it :-)
	*/
	// typedef matrix (std::vector< std::vector<double> >);

	__device__ void  cross
		(
		double vec1[3], double vec2[3], double outvec[3]
		);

	__device__ void  addvec
		(
		double a1, double vec1[3],
		double a2, double vec2[3],
		double vec3[3]
		);

	__device__ void  matvecmult
		(
		double mat[3][3],
		double vec[3],
		double vecout[3]
		);

	void  matmult
		(
		double** mat1,
		double** mat2,
		double** mat3,
		int mat1r, int mat1c, int mat2c
		);

	__device__ void  mattrans
		(
		double mat1[3][3],
		double mat2[3][3],
		int mat1r, int mat1c
		);

};  // namespace

#endif
