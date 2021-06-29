/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef CORE_MATRIX2D_H_
#define CORE_MATRIX2D_H_

#include <complex>
#include "xmipp_random_mode.h"
#include "xmipp_macros.h"
#include "xmipp_error.h"
#include "xmipp_strings.h"

class FileName;

#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif

template<typename T>
class Matrix1D;

// Forward declarations
template<typename T>
class Matrix2D;

template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d);

template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b);

template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v);

void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);


/** Cholesky decomposition.
 * Given M, this function decomposes M as M=L*L^t where L is a lower triangular matrix.
 * M must be positive semi-definite.
 */
void cholesky(const Matrix2D<double> &M, Matrix2D<double> &L);

/** Schur decomposition.
 * Given M, this function decomposes M as M = O*T*O' where O is an orthogonal matrix.
 */
void schur(const Matrix2D<double> &M, Matrix2D<double> &O, Matrix2D<double> &T);

/** @defgroup Matrices Matrix2D Matrices
 * @ingroup DataLibrary
 */
//@{
/** @name Matrices speed up macros */
//@{

/** Array access.
 *
 * This macro gives you access to the array (T)
 */
#define MATRIX2D_ARRAY(m) ((m).mdata)

/** For all elements in the array
 *
 * This macro is used to generate loops for the matrix in an easy way. It
 * defines internal indexes 'i' and 'j' which ranges the matrix using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
 * {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D(m) \
    for (size_t i=0; i<(m).mdimy; i++) \
        for (size_t j=0; j<(m).mdimx; j++)

/** Access to a matrix element
 * v is the array, i and j define the element v_ij.
 *
  * @code
 * MAT_ELEM(m, 0, 0) = 1;
 * val = MAT_ELEM(m, 0, 0);
 * @endcode
 */
#define MAT_ELEM(m,i,j) ((m).mdata[((size_t)i)*(m).mdimx+((size_t)j)])

/** X dimension of the matrix
 */
#define MAT_XSIZE(m) ((m).mdimx)

/** Y dimension of the matrix
 */
#define MAT_YSIZE(m) ((m).mdimy)

/** Total elements of the matrix
 */
#define MAT_SIZE(m) ((m).mdim)

/** Matrix element: Element access
 *
 * This is just a redefinition
 * of the function above
 * @code
 * dMij(m, -2, 1) = 1;
 * val = dMij(m, -2, 1);
 * @endcode
 */
#define dMij(m, i, j)  MAT_ELEM(m, i, j)

/** Matrix element: Element access
 *
 * This is just a redefinition
 * of the function above
 */
#define dMn(m, n)  ((m).mdata[(n)])

/** Matrix (3x3) by vector (3x1) (a=M*b)
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropriate size. You can reuse the vector b to store the results (that
 * is, M3x3_BY_V3x1(b, M, b);, is allowed).
 *
 * @code
 * double example
 * {
 *     SPEED_UP_temps;
 *
 *     Matrix1D< double > a(3), b(3);
 *     Matrix2D< double > M(3, 3);
 *
 *     M.init_random(0, 1);
 *     b.init_random(0, 1);
 *     M3x3_BY_V3x1(a, M, b);
 *
 *     return a.sum();
 * }
 * @endcode
 */
#define M3x3_BY_V3x1(a, M, b) { \
        spduptmp0 = dMn(M, 0) * XX(b) + dMn(M, 1) * YY(b) + dMn(M, 2) * ZZ(b); \
        spduptmp1 = dMn(M, 3) * XX(b) + dMn(M, 4) * YY(b) + dMn(M, 5) * ZZ(b); \
        spduptmp2 = dMn(M, 6) * XX(b) + dMn(M, 7) * YY(b) + dMn(M, 8) * ZZ(b); \
        XX(a) = spduptmp0; YY(a) = spduptmp1; ZZ(a) = spduptmp2; }

/** Matrix (3x3) by Matrix (3x3) (A=B*C)
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropriate size. You can reuse any of the multiplicands to store the
 * results (that is, M3x3_BY_M3x3(A, A, B);, is allowed).
 */
#define M3x3_BY_M3x3(A, B, C) { \
        spduptmp0 = dMn(B,0) * dMn(C,0) + dMn(B,1) * dMn(C,3) + dMn(B,2) * dMn(C,6); \
        spduptmp1 = dMn(B,0) * dMn(C,1) + dMn(B,1) * dMn(C,4) + dMn(B,2) * dMn(C,7); \
        spduptmp2 = dMn(B,0) * dMn(C,2) + dMn(B,1) * dMn(C,5) + dMn(B,2) * dMn(C,8); \
        spduptmp3 = dMn(B,3) * dMn(C,0) + dMn(B,4) * dMn(C,3) + dMn(B,5) * dMn(C,6); \
        spduptmp4 = dMn(B,3) * dMn(C,1) + dMn(B,4) * dMn(C,4) + dMn(B,5) * dMn(C,7); \
        spduptmp5 = dMn(B,3) * dMn(C,2) + dMn(B,4) * dMn(C,5) + dMn(B,5) * dMn(C,8); \
        spduptmp6 = dMn(B,6) * dMn(C,0) + dMn(B,7) * dMn(C,3) + dMn(B,8) * dMn(C,6); \
        spduptmp7 = dMn(B,6) * dMn(C,1) + dMn(B,7) * dMn(C,4) + dMn(B,8) * dMn(C,7); \
        spduptmp8 = dMn(B,6) * dMn(C,2) + dMn(B,7) * dMn(C,5) + dMn(B,8) * dMn(C,8); \
        dMn(A, 0) = spduptmp0; \
        dMn(A, 1) = spduptmp1; \
        dMn(A, 2) = spduptmp2; \
        dMn(A, 3) = spduptmp3; \
        dMn(A, 4) = spduptmp4; \
        dMn(A, 5) = spduptmp5; \
        dMn(A, 6) = spduptmp6; \
        dMn(A, 7) = spduptmp7; \
        dMn(A, 8) = spduptmp8; }

/** Matrix (2x2) by vector (2x1) (a=M*b)
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropriate size. You can reuse the vector b to store the results (that
 * is, M2x2_BY_V2x1(b, M, b);, is allowed).
 *
 * @code
 * double example
 * {
 *     SPEED_UP_temps;
 *
 *     Matrix1D< double > a(2), b(2);
 *     Matrix2D< double > M(2, 2);
 *
 *     M.init_random(0, 1);
 *     b.init_random(0, 1);
 *
 *     M2x2_BY_V2x1(a, M, b);
 *
 *     return a.sum();
 * }
 * @endcode
 */
#define M2x2_BY_V2x1(a, M, b) { \
        spduptmp0 = dMn(M, 0) * XX(b) + dMn(M, 1) * YY(b); \
        spduptmp1 = dMn(M, 2) * XX(b) + dMn(M, 3) * YY(b); \
        XX(a) = spduptmp0; \
        YY(a) = spduptmp1; }

/** Matrix (2x2) by constant (M2=M1*k)
 *
 * You must create the result matrix with the appropriate size. You can reuse
 * the matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is
 * allowed).
 */
#define M2x2_BY_CT(M2, M1, k) { \
        dMn(M2, 0) = dMn(M1, 0) * k; \
        dMn(M2, 1) = dMn(M1, 1) * k; \
        dMn(M2, 2) = dMn(M1, 2) * k; \
        dMn(M2, 3) = dMn(M1, 3) * k; }

/** Matrix (3x3) by constant (M2=M1*k)
 *
 * You must create the result matrix with the appropriate size. You can reuse the
 * matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is allowed).
 */
#define M3x3_BY_CT(M2, M1, k) { \
        dMn(M2, 0) = dMn(M1, 0) * k; \
        dMn(M2, 1) = dMn(M1, 1) * k; \
        dMn(M2, 2) = dMn(M1, 2) * k; \
        dMn(M2, 3) = dMn(M1, 3) * k; \
        dMn(M2, 4) = dMn(M1, 4) * k; \
        dMn(M2, 5) = dMn(M1, 5) * k; \
        dMn(M2, 6) = dMn(M1, 6) * k; \
        dMn(M2, 7) = dMn(M1, 7) * k; \
        dMn(M2, 8) = dMn(M1, 8) * k; }
/** Matrix (4x4) by constant (M2=M1*k)
 *
 * You must create the result matrix with the appropriate size. You can reuse the
 * matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is allowed).
 */
#define M4x4_BY_CT(M2, M1, k) { \
        dMn(M2, 0) = dMn(M1, 0) * k; \
        dMn(M2, 1) = dMn(M1, 1) * k; \
        dMn(M2, 2) = dMn(M1, 2) * k; \
        dMn(M2, 3) = dMn(M1, 3) * k; \
        dMn(M2, 4) = dMn(M1, 4) * k; \
        dMn(M2, 5) = dMn(M1, 5) * k; \
        dMn(M2, 6) = dMn(M1, 6) * k; \
        dMn(M2, 7) = dMn(M1, 7) * k; \
        dMn(M2, 8) = dMn(M1, 8) * k; \
  dMn(M2, 9) = dMn(M1, 9) * k; \
  dMn(M2,10) = dMn(M1,10) * k; \
  dMn(M2,11) = dMn(M1,11) * k; \
  dMn(M2,12) = dMn(M1,12) * k; \
  dMn(M2,13) = dMn(M1,13) * k; \
  dMn(M2,14) = dMn(M1,14) * k; \
  dMn(M2,15) = dMn(M1,15) * k;}

/** Inverse of a matrix (2x2)
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M2x2_INV(Ainv, A) { \
  spduptmp0 = dMn(A,0) * dMn(A,3) - dMn(A,1) * dMn(A,2);\
  if (spduptmp0==0.0) \
   REPORT_ERROR(ERR_NUMERICAL,"2x2 matrix is not invertible"); \
        spduptmp0 = 1.0 / spduptmp0; \
        dMn(Ainv, 0) =  dMn(A,3); \
        dMn(Ainv, 1) = -dMn(A,1); \
        dMn(Ainv, 2) = -dMn(A,2); \
        dMn(Ainv, 3) =  dMn(A,0); \
        M2x2_BY_CT(Ainv, Ainv, spduptmp0); }

/** Inverse of a matrix (3x3)
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M3x3_INV(Ainv, A) { \
        dMn(Ainv, 0) =   dMn(A,8)*dMn(A,4)-dMn(A,7)*dMn(A,5); \
        dMn(Ainv, 1) = -(dMn(A,8)*dMn(A,1)-dMn(A,7)*dMn(A,2)); \
        dMn(Ainv, 2) =   dMn(A,5)*dMn(A,1)-dMn(A,4)*dMn(A,2); \
        dMn(Ainv, 3) = -(dMn(A,8)*dMn(A,3)-dMn(A,6)*dMn(A,5)); \
        dMn(Ainv, 4) =   dMn(A,8)*dMn(A,0)-dMn(A,6)*dMn(A,2); \
        dMn(Ainv, 5) = -(dMn(A,5)*dMn(A,0)-dMn(A,3)*dMn(A,2)); \
        dMn(Ainv, 6) =   dMn(A,7)*dMn(A,3)-dMn(A,6)*dMn(A,4); \
        dMn(Ainv, 7) = -(dMn(A,7)*dMn(A,0)-dMn(A,6)*dMn(A,1)); \
        dMn(Ainv, 8) =   dMn(A,4)*dMn(A,0)-dMn(A,3)*dMn(A,1); \
        spduptmp0 = dMn(A,0)*dMn(Ainv,0)+dMn(A,3)*dMn(Ainv,1)+dMn(A,6)*dMn(Ainv,2); \
  if (spduptmp0==0.0) \
   REPORT_ERROR(ERR_NUMERICAL,"3x3 matrix is not invertible"); \
     spduptmp0 = 1.0 / spduptmp0; \
        M3x3_BY_CT(Ainv, Ainv, spduptmp0); }

/** Inverse of a matrix (4x4)
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M4x4_INV(Ainv, A) { \
        dMn(Ainv, 0) =   dMn(A,5)*(dMn(A,10)*dMn(A,15)-dMn(A,11)*dMn(A,14))+\
             dMn(A,6)*(dMn(A,11)*dMn(A,13)-dMn(A,9) *dMn(A,15))+\
             dMn(A,7)*(dMn(A,9) *dMn(A,14)-dMn(A,10)*dMn(A,13));\
  dMn(Ainv, 1) =   dMn(A,1)*(dMn(A,11)*dMn(A,14)-dMn(A,10)*dMn(A,15))+\
       dMn(A,2)*(dMn(A,9) *dMn(A,15)-dMn(A,11)*dMn(A,13))+\
       dMn(A,3)*(dMn(A,10)*dMn(A,13)-dMn(A,9) *dMn(A,14));\
  dMn(Ainv, 2) =   dMn(A,1)*(dMn(A,6) *dMn(A,15)-dMn(A,7) *dMn(A,14))+\
       dMn(A,2)*(dMn(A,7) *dMn(A,13)-dMn(A,5) *dMn(A,15))+\
       dMn(A,3)*(dMn(A,5) *dMn(A,14)-dMn(A,6) *dMn(A,13));\
  dMn(Ainv, 3) =   dMn(A,1)*(dMn(A,7) *dMn(A,10)-dMn(A,6) *dMn(A,11))+\
       dMn(A,2)*(dMn(A,5) *dMn(A,11)-dMn(A,7) *dMn(A,9))+\
       dMn(A,3)*(dMn(A,6) *dMn(A,9) -dMn(A,5) *dMn(A,10));\
  dMn(Ainv, 4) =   dMn(A,4)*(dMn(A,11)*dMn(A,14)-dMn(A,10)*dMn(A,15))+\
       dMn(A,6)*(dMn(A,8) *dMn(A,15)-dMn(A,11)*dMn(A,12))+\
       dMn(A,7)*(dMn(A,10)*dMn(A,12)-dMn(A,8) *dMn(A,14));\
  dMn(Ainv, 5) =   dMn(A,0)*(dMn(A,10)*dMn(A,15)-dMn(A,11)*dMn(A,14))+\
       dMn(A,2)*(dMn(A,11)*dMn(A,12)-dMn(A,8) *dMn(A,15))+\
       dMn(A,3)*(dMn(A,8) *dMn(A,14)-dMn(A,10)*dMn(A,12));\
  dMn(Ainv, 6) =   dMn(A,0)*(dMn(A,7) *dMn(A,14)-dMn(A,6) *dMn(A,15))+\
       dMn(A,2)*(dMn(A,4) *dMn(A,15)-dMn(A,7) *dMn(A,12))+\
       dMn(A,3)*(dMn(A,6) *dMn(A,12)-dMn(A,4) *dMn(A,14));\
  dMn(Ainv, 7) =   dMn(A,0)*(dMn(A,6) *dMn(A,11)-dMn(A,7) *dMn(A,10))+\
       dMn(A,2)*(dMn(A,7) *dMn(A,8) -dMn(A,4) *dMn(A,11))+\
       dMn(A,3)*(dMn(A,4) *dMn(A,10)-dMn(A,6) *dMn(A,8));\
  dMn(Ainv, 8) =   dMn(A,4)*(dMn(A,9) *dMn(A,15)-dMn(A,11)*dMn(A,13))+\
       dMn(A,5)*(dMn(A,11)*dMn(A,12)-dMn(A,8) *dMn(A,15))+\
       dMn(A,7)*(dMn(A,8) *dMn(A,13)-dMn(A,9) *dMn(A,12));\
  dMn(Ainv, 9) =   dMn(A,0)*(dMn(A,11)*dMn(A,13)-dMn(A,9) *dMn(A,15))+\
       dMn(A,1)*(dMn(A,8) *dMn(A,15)-dMn(A,11)*dMn(A,12))+\
       dMn(A,3)*(dMn(A,9) *dMn(A,12)-dMn(A,8) *dMn(A,13));\
  dMn(Ainv,10) =   dMn(A,0)*(dMn(A,5) *dMn(A,15)-dMn(A,7) *dMn(A,13))+\
       dMn(A,1)*(dMn(A,7) *dMn(A,12)-dMn(A,4) *dMn(A,15))+\
       dMn(A,3)*(dMn(A,4) *dMn(A,13)-dMn(A,5) *dMn(A,12));\
  dMn(Ainv,11) =   dMn(A,0)*(dMn(A,7) *dMn(A,9) -dMn(A,5) *dMn(A,11))+\
       dMn(A,1)*(dMn(A,4) *dMn(A,11)-dMn(A,7) *dMn(A,8))+\
       dMn(A,3)*(dMn(A,5) *dMn(A,8) -dMn(A,4) *dMn(A,9));\
  dMn(Ainv,12) =   dMn(A,4)*(dMn(A,10)*dMn(A,13)-dMn(A,9) *dMn(A,14))+\
       dMn(A,5)*(dMn(A,8) *dMn(A,14)-dMn(A,10)*dMn(A,12))+\
       dMn(A,6)*(dMn(A,9) *dMn(A,12)-dMn(A,8) *dMn(A,13));\
  dMn(Ainv,13) =   dMn(A,0)*(dMn(A,9) *dMn(A,14)-dMn(A,10)*dMn(A,13))+\
       dMn(A,1)*(dMn(A,10)*dMn(A,12)-dMn(A,8) *dMn(A,14))+\
       dMn(A,2)*(dMn(A,8) *dMn(A,13)-dMn(A,9) *dMn(A,12));\
  dMn(Ainv,14) =   dMn(A,0)*(dMn(A,6) *dMn(A,13)-dMn(A,5) *dMn(A,14))+\
       dMn(A,1)*(dMn(A,4) *dMn(A,14)-dMn(A,6) *dMn(A,12))+\
       dMn(A,2)*(dMn(A,5) *dMn(A,12)-dMn(A,4) *dMn(A,13));\
  dMn(Ainv,15) =   dMn(A,0)*(dMn(A,5) *dMn(A,10)-dMn(A,6) *dMn(A,9))+\
       dMn(A,1)*(dMn(A,6) *dMn(A,8) -dMn(A,4) *dMn(A,10))+\
       dMn(A,2)*(dMn(A,4) *dMn(A,9) -dMn(A,5) *dMn(A,8));\
        spduptmp0 = dMn(A,0)*(dMn(A,5)*(dMn(A,10)*dMn(A,15)-dMn(A,11)*dMn(A,14))\
                              +dMn(A,6)*(dMn(A,11)*dMn(A,13)-dMn(A,9) *dMn(A,15))\
                              +dMn(A,7)*(dMn(A,9) *dMn(A,14)-dMn(A,10)*dMn(A,13)))\
                    +dMn(A,1)*(dMn(A,4)*(dMn(A,11)*dMn(A,14)-dMn(A,10)*dMn(A,15))\
                              +dMn(A,6)*(dMn(A,8) *dMn(A,15)-dMn(A,11)*dMn(A,12))\
                              +dMn(A,7)*(dMn(A,10)*dMn(A,12)-dMn(A,8) *dMn(A,14)))\
        +dMn(A,2)*(dMn(A,4)*(dMn(A,9) *dMn(A,15)-dMn(A,11)*dMn(A,13))\
         +dMn(A,5)*(dMn(A,11)*dMn(A,12)-dMn(A,8) *dMn(A,15))\
         +dMn(A,7)*(dMn(A,8) *dMn(A,13)-dMn(A,9) *dMn(A,12)))\
        +dMn(A,3)*(dMn(A,4)*(dMn(A,10)*dMn(A,13)-dMn(A,9) *dMn(A,14))\
         +dMn(A,5)*(dMn(A,8) *dMn(A,14)-dMn(A,10)*dMn(A,12))\
         +dMn(A,6)*(dMn(A,9) *dMn(A,12)-dMn(A,8) *dMn(A,13))); \
  if (spduptmp0==0.0) \
   REPORT_ERROR(ERR_NUMERICAL,"4x4 matrix is not invertible"); \
  spduptmp0 = 1.0 / spduptmp0; \
        M4x4_BY_CT(Ainv, Ainv, spduptmp0); }

/** Matrix2D class */
template<typename T>
class Matrix2D
{
public:
    // The array itself
    T* mdata;

    // Destroy data
    bool destroyData;

    // Mapped data
    bool mappedData;

    // File descriptor for mapped files
    int fdMap;

    // Mapped data original pointer
    char* mdataOriginal;

    // Number of elements in X
    size_t mdimx;

    // Number of elements in Y
    size_t mdimy;

    // Total number of elements
    size_t mdim;
    //@}

    /// @name Constructors
    //@{
    /** Empty constructor
     */
    Matrix2D()
    {
        coreInit();
    }

    Matrix2D(const FileName &fnMappedMatrix, int Ydim, int Xdim, size_t offset=0)
    {
        coreInit(fnMappedMatrix,Ydim,Xdim,offset);
    }

    /** Dimension constructor
     */
    Matrix2D(int Ydim, int Xdim)
    {
        coreInit();
        initZeros(Ydim, Xdim);
    }

    /** Copy constructor
     */
    Matrix2D(const Matrix2D<T>& v)
    {
        coreInit();
        *this = v;
    }

    /** Destructor.
     */
    ~Matrix2D()
    {
        coreDeallocate();
    }

    /** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix2D<T>& operator=(const Matrix2D<T>& op1);
    //@}

    /// @name Core memory operations for Matrix2D
    //@{
    /** Clear.
     */
    void clear()
    {
        coreDeallocate();
        coreInit();
    }

    /** Method will convert Matrix2D matrix to float[3][3] */
    inline void convertTo(float out[3][3]) const {
    	for (int i = 0; i < 3; i++) {
    		for (int j = 0; j < 3; j++) {
    			out[i][j] = (*this)(i, j);
    		}
    	}
    }

    /** Core init from mapped file.
     * Offset is in bytes. */
    void coreInit(const FileName &fn, int Ydim, int Xdim, size_t offset=0);

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit();

    /** Core allocate.
     */
    void coreAllocate( int _mdimy, int _mdimx);

    /** Core deallocate.
     * Free all mdata.
     */
    void coreDeallocate();
    //@}

    /// @name Size and shape of Matrix2D
    //@{
    /** Resize to a given size
     */
    void resize(size_t Ydim, size_t Xdim, bool noCopy=false);

    /** Resize according to a pattern.
     *
     * This function resize the actual array to the same size and origin
     * as the input pattern. If the actual array is larger than the pattern
     * then the trailing values are lost, if it is smaller then 0's are
     * added at the end
     *
     * @code
     * v2.resize(v1);
     * // v2 has got now the same structure as v1
     * @endcode
     */
    template<typename T1>
    void resize(const Matrix2D<T1> &v)
    {
        if (mdimx != v.mdimx || mdimy != v.mdimy)
            resize(v.mdimy, v.mdimx);
    }

    /** Resize to a given size (don't copy old elements)
     */
    inline void resizeNoCopy(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim, true);
    }

    /** Resize according to a pattern.
     *  Do not copy old elements.
     */
    template<typename T1>
    inline void resizeNoCopy(const Matrix2D<T1> &v)
    {
        if (mdimx != v.mdimx || mdimy != v.mdimy)
            resize(v.mdimy, v.mdimx, true);
    }

    /** Map to file.
     * The matrix is mapped to a file. The file is presumed to be already created with
     * enough space for a matrix of size Ydim x Xdim.
     * Offset is in bytes.
     */
    void mapToFile(const FileName &fn, int Ydim, int Xdim, size_t offset=0);

    /** Extract submatrix and assign to this object.
     */
    void submatrix(int i0, int j0, int iF, int jF);

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    inline bool sameShape(const Matrix2D<T1>& op) const
    {
        return ((mdimx == op.mdimx) && (mdimy == op.mdimy));
    }

    /** X dimension
     *
     * Returns X dimension
     */
    inline size_t Xdim() const
    {
        return mdimx;
    }

    /** Y dimension
     *
     * Returns Y dimension
     */
    inline size_t Ydim() const
    {
        return mdimy;
    }
    //@}

    /// @name Initialization of Matrix2D values
    //@{
    /** Same value in all components.
     *
     * The constant must be of a type compatible with the array type, ie,
     * you cannot  assign a double to an integer array without a casting.
     * It is not an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initConstant(3.14);
     * @endcode
     */
    inline void initConstant(T val)
    {
        for (size_t j = 0; j < mdim; j++)
            mdata[j] = val;
    }

    /** Initialize to zeros with a given size.
     */
    inline void initConstant(size_t Ydim, size_t Xdim, T val)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        initConstant(val);
    }

    /** Initialize to zeros with current size.
     *
     * All values are set to 0. The current size and origin are kept. It is not
     * an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initZeros();
     * @endcode
     */
    inline void initZeros()
    {
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to zeros with a given size.
     */
    inline void initZeros(size_t Ydim, size_t Xdim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to zeros following a pattern.
      *
      * All values are set to 0, and the origin and size of the pattern are
      * adopted.
      *
      * @code
      * v2.initZeros(v1);
      * @endcode
      */
    template <typename T1>
    void initZeros(const Matrix2D<T1>& op)
    {
        if (mdimx!=op.mdimx || mdimy!=op.mdimy)
            resizeNoCopy(op);
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to random random numbers, uniform or gaussian
     */
    void initRandom(size_t Ydim, size_t Xdim, double op1, double op2, RandomMode mode = RND_UNIFORM);

    /** Initialize to gaussian numbers */
    void initGaussian(int Ydim, int Xdim, double op1=0., double op2=1.);

    /** 2D Identity matrix of current size
     *
     * If actually the matrix is not squared then an identity matrix is
     * generated of size (Xdim x Xdim).
     *
     * @code
     * m.initIdentity();
     * @endcode
     */
    inline void initIdentity()
    {
        initIdentity(MAT_XSIZE(*this));
    }

    /** 2D Identity matrix of a given size
     *
     * A (dim x dim) identity matrix is generated.
     *
     * @code
     * m.initIdentity(3);
     * @endcode
     */
    inline void initIdentity(int dim)
    {
        initZeros(dim, dim);
        for (int i = 0; i < dim; i++)
            MAT_ELEM(*this,i,i) = 1;
    }
    /** 2D gaussian matrix of a given size and with a given variance. The amplitude of the Gaussian is set to 1.
     *
     * A (dim x dim) gaussian matrix is generated.
     *
     * @code
     * m.initGaussian(3,1);
     * @endcode
     */
    void initGaussian(int dim, double var);
    //@}

    /// @name Operators for Matrix2D
    //@{

    /** Matrix element access
     */
    inline T& operator()(int i, int j) const
    {
        return MAT_ELEM((*this),i,j);
    }
    /** Parenthesis operator for phyton
    */
    inline void setVal(T val,int y, int x)
    {
        MAT_ELEM((*this),y,x)=val;
    }
    /** Parenthesis operator for phyton
    */
    inline T getVal( int y, int x) const
    {
        return MAT_ELEM((*this),y,x);
    }

    /** v3 = v1 * k.
     */
    inline Matrix2D<T> operator*(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (size_t i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] * op1;
        return tmp;
    }

    /** v3 = v1 / k.
     */
    inline Matrix2D<T> operator/(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (size_t i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] / op1;
        return tmp;
    }

    /** v3 = k * v2.
     */
    inline friend Matrix2D<T> operator*(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp(op2);
        for (size_t i=0; i < op2.mdim; i++)
            tmp.mdata[i] = op1 * op2.mdata[i];
        return tmp;
    }

    /** v3 *= k.
      */
    inline void operator*=(T op1)
    {
        for (size_t i=0; i < mdim; i++)
            mdata[i] *= op1;
    }

    /** v3 /= k.
      */
    inline void operator/=(T op1)
    {
        for (size_t i=0; i < mdim; i++)
            mdata[i] /= op1;
    }

    /** Matrix by vector multiplication
     *
     * @code
     * v2 = A*v1;
     * @endcode
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const;

    /** Matrix by Matrix multiplication
     *
     * @code
     * C = A*B;
     * @endcode
     */
    Matrix2D<T> operator*(const Matrix2D<T>& op1) const;

    /** Matrix summation
     *
     * @code
     * C = A + B;
     * @endcode
     */
    Matrix2D<T> operator+(const Matrix2D<T>& op1) const;

    /** Matrix summation
     *
     * @code
     * A += B;
     * @endcode
     */
    void operator+=(const Matrix2D<T>& op1) const;

    /** Matrix subtraction
     *
     * @code
     * C = A - B;
     * @endcode
     */
    Matrix2D<T> operator-(const Matrix2D<T>& op1) const;

    /** Matrix subtraction
     *
     * @code
     * A -= B;
     * @endcode
     */
    void operator-=(const Matrix2D<T>& op1) const;
    /** Equality.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const Matrix2D<T>& op,
               double accuracy = XMIPP_EQUAL_ACCURACY) const;
    /** Equality.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy and without SIGN).
     */
    bool equalAbs(const Matrix2D<T>& op,
               double accuracy = XMIPP_EQUAL_ACCURACY) const;
    //@}

    /// @name Utilities for Matrix2D
    //@{
    /** Maximum of the values in the array.
      *
      * The returned value is of the same type as the type of the array.
      */
    T computeMax() const;

    /** Minimum of the values in the array.
       *
       * The returned value is of the same type as the type of the array.
       */
    T computeMin() const;

    /** Maximum and minimum of the values in the array. */
    void computeMaxAndMin(T &maxValue, T &minValue) const;

    /** Get row sum. */
    void rowSum(Matrix1D<T> &sum) const;

    /** Get column sum. */
        void colSum(Matrix1D<T> &sum) const;

    /** Get row energy sum.
     * Sum of the squared values by row */
    void rowEnergySum(Matrix1D<T> &sum) const;

    /** Produce a 2D array suitable for working with Numerical Recipes
    *
    * This function must be used only as a preparation for routines which need
    * that the first physical index is 1 and not 0 as it usually is in C. New
    * memory is needed to hold the new double pointer array.
    */
    T** adaptForNumericalRecipes() const;

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    inline T* adaptForNumericalRecipes2() const
    {
        return mdata - 1 - mdimx;
    }

    /** Load 2D array from numerical recipes result.
     */
    void loadFromNumericalRecipes(T** m, int Ydim, int Xdim);

    /** Kill a 2D array produced for numerical recipes
     *
     * The allocated memory is freed.
     */
    inline void killAdaptationForNumericalRecipes(T** m) const
    {
        free_Tmatrix(m, 1, mdimy, 1, mdimx);
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes2(T** m) const
        {}

    /** Read this matrix from file.
     *  The matrix is assumed to be already resized.
      */
    void read(const FileName &fn);

    /** Write this matrix to file
      */
    void write(const FileName &fn) const;
    /** Show matrix
      */
    friend std::ostream& operator<<(std::ostream& ostrm, const Matrix2D<T>& v)
    {
        if (v.Xdim() == 0 || v.Ydim() == 0)
            ostrm << "NULL matrix\n";
        else
        {
            ostrm << std::endl;
            double max_val = v.computeMax();
            int prec = bestPrecision(max_val, 10);

            for (size_t i = 0; i < v.Ydim(); i++)
            {
                for (size_t j = 0; j < v.Xdim(); j++)
                {
                    ostrm << floatToString((double) v(i, j), 10, prec) << ' ';
                }
                ostrm << std::endl;
            }
        }

        return ostrm;
    }

    /** Makes a matrix from a vector
     *
     * The origin of the matrix is set such that it has one of the index origins
     * (X or Y) to the same value as the vector, and the other set to 0
     * according to the shape.
     *
     * @code
     * Matrix2D< double > m = fromVector(v);
     * @endcode
     */
    void fromVector(const Matrix1D<T>& op1);

    /** Makes a vector from a matrix
     *
     * An exception is thrown if the matrix is not a single row or a single
     * column. The origin of the vector is set according to the one of the
     * matrix.
     *
     * @code
     * Matrix1D< double > v;
     * m.toVector(v);
     * @endcode
     */
    void toVector(Matrix1D<T>& op1) const;

    /**Copy matrix to stl::vector
     */
    void copyToVector(std::vector<T> &v)
    {
        v.assign(mdata, mdata+mdim);
    }
    /**Copy stl::vector to matrix
      */
    inline void copyFromVector(std::vector<T> &v,int Xdim, int Ydim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resizeNoCopy(Ydim, Xdim);
        copy( v.begin(), v.begin()+v.size(), mdata);
    }

    /** Get row
     *
     * This function returns a row vector corresponding to the chosen
     * row inside the nth 2D matrix, the numbering of the rows is also
     * logical not physical.
     *
     * @code
     * std::vector< double > v;
     * m.getRow(-2, v);
     * @endcode
     */
    void getRow(size_t i, Matrix1D<T>& v) const;

    /** Get Column
     *
     * This function returns a column vector corresponding to the
     * chosen column.
     *
     * @code
     * std::vector< double > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(size_t j, Matrix1D<T>& v) const;

    /** Set Row
     *
     * This function sets a row vector corresponding to the chosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(size_t i, const Matrix1D<T>& v);

    /** Set Column
     *
     * This function sets a column vector corresponding to the chosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
     * @endcode
     */
    void setCol(size_t j, const Matrix1D<T>& v);

    /** Compute row means */
    void computeRowMeans(Matrix1D<double> &Xmr) const;

    /** Compute row means */
    void computeColMeans(Matrix1D<double> &Xmr) const;

     /** Set constant column.
     * Set a given column to a constant value A(i,j)=val;
     */
    inline void setConstantCol(size_t j, T v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(ERR_MATRIX_EMPTY, "setCol: Target matrix is empty");

        if (j>= mdimx)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setCol: Matrix subscript (j) out of range");

        for (size_t i = 0; i < mdimy; i++)
            MAT_ELEM(*this,i, j) = v;
    }

    /** Get diagonal.
     * It is assumed that the matrix is squared
     */
    void getDiagonal(Matrix1D<T> &d) const;

    /** Trace
     * Sum of the values in the diagonal
     */
    T trace() const;

    /** Determinant of a matrix
     *
     * An exception is thrown if the matrix is not squared or it is empty.
     *
     * @code
     * double det = m.det();
     * @endcode
     */
    T det() const;

    ///determinat of 3x3 matrix
    inline T det3x3() const
    {
        return (
                   dMij(*this,0,0)*( dMij(*this,2,2)*dMij(*this,1,1)-
                                     dMij(*this,2,1)*dMij(*this,1,2) )-
                   dMij(*this,1,0)*( dMij(*this,2,2)*dMij(*this,0,1)-
                                     dMij(*this,2,1)*dMij(*this,0,2) )+
                   dMij(*this,2,0)*( dMij(*this,1,2)*dMij(*this,0,1)-
                                     dMij(*this,1,1)*dMij(*this,0,2) )
               );
    }

    /** Frobenius norm of a matrix */
    inline double norm()
    {
    	double sum=0.;
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
    	{
    		T aux=MAT_ELEM(*this,i,j);
    		sum+=aux*aux;
    	}
    	return sqrt(sum);
    }

    /** Algebraic transpose of a Matrix
     *
     * You can use the transpose in as complex expressions as you like. The
     * origin of the vector is not changed.
     *
     * @code
     * v2 = v1.transpose();
     * @endcode
     */
    Matrix2D<T> transpose() const;

    /** Inverse of a matrix
     *
     * The matrix is inverted using a SVD decomposition. In fact the
     * pseudoinverse is returned.
     *
     * @code
     * Matrix2D< double > m1_inv;
     * m1.inv(m1_inv);
     * @endcode
     */
    void inv(Matrix2D<T>& result) const;

    /** Inverse of a matrix
     *
     * The matrix is inverted using a AlgLib.
     * Set LU to use LU decomposition
     *
     * @code
     * Matrix2D< double > m1_inv;
     * m1.inv(m1_inv);
     * @endcode
     */
    void invAlgLib(Matrix2D<T>& result, bool use_lu=false) const;

    /** Perform SVD decomposition
     *  *this = U * W * V^t
     * */

    void svd(Matrix2D<double> &U, Matrix1D<double> &W, Matrix2D<double> &V) const
    {
      svdcmp(*this, U, W, V);
    }

    /** Perform SVD decomposition and add an index vector
     * with the descending order of singular values
     */
    void eigs(Matrix2D<double> &U, Matrix1D<double> &W, Matrix2D<double> &V, Matrix1D<int> &indexes) const;

    /** Inverse of a matrix
     */
    Matrix2D<T> inv() const
    {
        Matrix2D<T> result;
        inv(result);

        return result;
    }

    /** Inverse the current matrix
     */
    void selfInverse()
    {
        Matrix2D<T> auxMatrix(*this);
        auxMatrix.inv(*this);
    }

    /** True if the matrix is identity
     *
     * @code
     * if (m.isIdentity())
     *     std::cout << "The matrix is identity\n";
     * @endcode
     */
    bool isIdentity() const;
    //@}
};

typedef Matrix2D<double> DMatrix;
typedef Matrix2D<int> IMatrix;

template<typename T>
bool operator==(const Matrix2D<T>& op1, const Matrix2D<T>& op2)
{
    return op1.equal(op2);
}

/**@name Matrix Related functions
 * These functions are not methods of Matrix2D
 */
//@{
/** LU Decomposition
 */
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d);


/** LU Backsubstitution
 */
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b);

/** SVD Backsubstitution
 */
void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);

//#define VIA_NR
#define VIA_BILIB
/** SVD Decomposition
 */
template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v);

/** Generalized eigenvector decomposition.
 * Solves the problem Av=dBv.
 * The decomposition is such that A=B P D P^-1. A and B must be square matrices of the same size.
 */
void generalizedEigs(const Matrix2D<double> &A, const Matrix2D<double> &B, Matrix1D<double> &D, Matrix2D<double> &P);

/** First eigenvectors of a real, symmetric matrix.
 * Solves the problem Av=dv.
 * Only the eigenvectors of the largest M eigenvalues are returned as columns of P
 */
void firstEigs(const Matrix2D<double> &A, size_t M, Matrix1D<double> &D, Matrix2D<double> &P, bool Pneeded=true);

/** Last eigenvectors of a real, symmetric matrix.
 * Solves the problem Av=dv.
 * Only the eigenvectors of the smallest M eigenvalues are returned as columns of P
 */
void lastEigs(const Matrix2D<double> &A, size_t M, Matrix1D<double> &D, Matrix2D<double> &P);

/** Compute eigenvectors between two indexes of a real, symmetric matrix.
 * Solves the problem Av=dv.
 * Only the eigenvectors of the smallest eigenvalues between indexes I1 and I2 are returned as columns of P. Indexes start at 0.
 */
void eigsBetween(const Matrix2D<double> &A, size_t I1, size_t I2, Matrix1D<double> &D, Matrix2D<double> &P);

/** Compute all eigenvalues, even if they are complex */
void allEigs(const Matrix2D<double> &A, std::vector< std::complex<double> > &eigs);

/** Find connected components of a graph.
 * Assuming that the matrix G represents an undirected graph (of size NxN), this function returns a vector (of size N) that indicates
 * for each element which is the number of its connected component.
 */
void connectedComponentsOfUndirectedGraph(const Matrix2D<double> &G, Matrix1D<int> &component);

/** Conversion from one type to another.
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T2>& v2)
{
    if (v1.mdim == 0)
    {
        v2.clear();
        return;
    }

    if (v1.mdimx!=v2.mdimx || v1.mdimy!=v2.mdimy)
        v2.resizeNoCopy(v1);
    for (unsigned long int n = 0; n < v1.mdim; n++)
        v2.mdata[n] = static_cast< T2 > (v1.mdata[n]);
}

/** Conversion from one type to another.
 * In some cases, the two types are the same. So a faster way is simply by assignment.
 */
template<typename T1>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T1>& v2)
{
    v2=v1;
}

/** Gram Schmidt orthogonalization by columns */
void orthogonalizeColumnsGramSchmidt(Matrix2D<double> &M);

/** Normalize columns.
 * So that they have zero mean and unit variance.
 */
void normalizeColumns(Matrix2D<double> &A);

/** Normalize columns.
 * So that the minimum is 0 and the maximum is 1
 */
void normalizeColumnsBetween0and1(Matrix2D<double> &A);

/** Subtract mean of columns.
 * So that they have zero mean.
 */
void subtractColumnMeans(Matrix2D<double> &A);

/** Matrix operation: B=A^t*A. */
void matrixOperation_AtA(const Matrix2D <double> &A, Matrix2D<double> &B);

/** Matrix operation: C=A*A^t. */
void matrixOperation_AAt(const Matrix2D <double> &A, Matrix2D<double> &C);

/** Matrix operation: C=A*B. */
void matrixOperation_AB(const Matrix2D <double> &A, const Matrix2D<double> &B, Matrix2D<double> &C);

/** Matrix operation: y=A*x. */
void matrixOperation_Ax(const Matrix2D <double> &A, const Matrix1D<double> &x, Matrix1D<double> &y);

/** Matrix operation: C=A*B^t. */
void matrixOperation_ABt(const Matrix2D <double> &A, const Matrix2D <double> &B, Matrix2D<double> &C);

/** Matrix operation: C=A^t*B. */
void matrixOperation_AtB(const Matrix2D <double> &A, const Matrix2D<double> &B, Matrix2D<double> &C);

/** Matrix operation: y=A^t*x. */
void matrixOperation_Atx(const Matrix2D <double> &A, const Matrix1D<double> &x, Matrix1D<double> &y);

/** Matrix operation: C=A^t*Bt. */
void matrixOperation_AtBt(const Matrix2D <double> &A, const Matrix2D<double> &B, Matrix2D<double> &C);

/** Matrix operation: B=X^t*A*X.
 * We know that the result B must be symmetric */
void matrixOperation_XtAX_symmetric(const Matrix2D<double> &X, const Matrix2D<double> &A, Matrix2D<double> &B);

/** Matrix operation: A=I+A */
void matrixOperation_IplusA(Matrix2D<double> &A);

/** Matrix operation: A=I-A */
void matrixOperation_IminusA(Matrix2D<double> &A);

/** Erase first column */
void eraseFirstColumn(Matrix2D<double> &A);

/** Keep columns between j0 and jF */
void keepColumns(Matrix2D<double> &A, int j0, int jF);

//@}
//@}

#endif /* MATRIX2D_H_ */
