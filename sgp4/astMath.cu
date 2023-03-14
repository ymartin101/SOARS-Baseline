/*     ----------------------------------------------------------------
*
*                                astMath.cu
*
*    this file contains miscellaneous math functions.
*
*                          Companion code for
*             Fundamentals of Astrodynamics and Applications
*                                  2013
*                            by David Vallado
*
*       (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
*
*    current :
*              11 jan 18  david vallado
*                           misc cleanup
*    changes :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               7 may 08  david vallado
*                           misc updates, fix sgn, show both options for matrix
*                           multiplications
*              22 jan 08  david vallado
*                           fix some minor errors, fixes to matmult
*              19 oct 07  david vallado
*                           fix sgn baseline comparison
*              30 may 07  david vallado
*                           3rd edition baseline
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */

#include "astMath.cuh"

namespace astMath
{

	/* -----------------------------------------------------------------------------
	*
	*                           procedure addvec
	*
	*  this procedure adds two vectors possibly multiplied by a constant.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    a1          - constant multiplier
	*    a2          - constant multiplier
	*    vec1        - vector number 1
	*    vec2        - vector number 2
	*
	*  outputs       :
	*    outvec      - vector result of a + b
	*
	*  locals        :
	*    row         - index
	*
	*  coupling      :
	*     none
	* --------------------------------------------------------------------------- */

		__device__ void    addvec
			(
			double a1, double vec1[3],
			double a2, double vec2[3],
			double vec3[3]
			)
		{
			int row;

			for (row = 0; row <= 2; row++)
			{
				vec3[row] = 0.0;
				vec3[row] = a1* vec1[row] + a2* vec2[row];
			}
		}  // addvec

/* -----------------------------------------------------------------------------
*
*                           procedure cross
*
*  this procedure crosses two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    outvec      - vector result of a x b
*
*  locals        :
*    none.
*
*  coupling      :
*    none
 ---------------------------------------------------------------------------- */

	__device__ void    cross
		(
		double vec1[3], double vec2[3], double outvec[3]
		)
	{
		outvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
		outvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
		outvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
	}  // cross

/* -----------------------------------------------------------------------------
*
*                           procedure matvecmult
*
*  this procedure multiplies a 3x3 matrix and a 3x1 vector together.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    mat         - 3 x 3 matrix
*    vec         - vector
*
*  outputs       :
*    vecout      - vector result of mat * vec
*
*  locals        :
*    row         - row index
*    col         - column index
*    ktr         - index
*
*  coupling      :
* --------------------------------------------------------------------------- */

	__device__ void    matvecmult
		(
		double mat[3][3],
		//          double mat[3][3],
		double vec[3],
		double vecout[3]
		)
	{
		int row, ktr;

		for (row = 0; row <= 2; row++)
		{
			vecout[row] = 0.0;
			for (ktr = 0; ktr <= 2; ktr++)
				vecout[row] = vecout[row] + mat[row][ktr] * vec[ktr];
		}
	}  // matvecmult

/* -----------------------------------------------------------------------------
*
*                           procedure matmult
*
*  this procedure multiplies two matricies up to 10x10 together.
*
*  author        : david vallado                  719-573-2600    7 dec 2007
*
*  inputs          description                    range / units
*    mat1        - matrix number 1
*    mat2        - matrix number 2
*    mat1r       - matrix number 1 rows
*    mat1c       - matrix number 1 columns
*    mat2c       - matrix number 2 columns
*
*  outputs       :
*    mat3        - matrix result of mat1 * mat2 of size mat1r x mat2c
*
*  locals        :
*    row         - row index
*    col         - column index
*    ktr         - index
*
*  coupling      :
* --------------------------------------------------------------------------- */

	void    matmult
		(
		double** mat1,
		double** mat2,
		double** mat3,
		//          double mat1[3][3],
		//          double mat2[3][3],
		//          double mat3[3][3],
		int mat1r, int mat1c, int mat2c
		)
	{
		int row, col, ktr;

		for (row = 0; row < mat1r; row++)
		{
			for (col = 0; col < mat2c; col++)
			{
				mat3[row][col] = 0.0;
				for (ktr = 0; ktr < mat1c; ktr++)
					mat3[row][col] = mat3[row][col] + mat1[row][ktr] * mat2[ktr][col];
			}
		}
	}  // matmult

/* -----------------------------------------------------------------------------
*
*                           procedure mattrans
*
*  this procedure finds the transpose of a matrix.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    mat1        - matrix number 1
*    mat1r       - matrix number 1 rows
*    mat1c       - matrix number 1 columns
*
*  outputs       :
*    mat2        - matrix result of transpose mat2
*
*  locals        :
*    row         - row index
*    col         - column index
*
*  coupling      :
* --------------------------------------------------------------------------- */

	__device__ void  mattrans
		(
		double mat1[3][3],
		double mat2[3][3],
		int mat1r, int mat1c
		)
	{
		int row, col;

		for (row = 0; row < mat1r; row++)
		{
			for (col = 0; col < mat1c; col++)
				mat2[col][row] = mat1[row][col];
		}
	}  // mattrans

}  // namespace
