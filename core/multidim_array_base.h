/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
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

#ifndef XMIPPCORE_CORE_MULTIDIM_ARRAY_BASE_H_
#define XMIPPCORE_CORE_MULTIDIM_ARRAY_BASE_H_

#include <stddef.h>
#include <iostream>
#include "xmipp_array_dim.h"
#include "xmipp_array_coord.h"

template<typename T>
class Matrix1D;

/** @name MultidimArraysSpeedUp Speed up macros
 *
 * This macros are defined to allow high speed in critical parts of your
 * program. They shouldn't be used systematically as usually there is no
 * checking on the correctness of the operation you are performing. Speed comes
 * from three facts: first, they are macros and no function call is performed
 * (although most of the critical functions are inline functions), there is no
 * checking on the correctness of the operation (it could be wrong and you are
 * not warned of it), and destination vectors are not returned saving time in
 * the copy constructor and in the creation/destruction of temporary vectors.
 */
//@{
/** Returns the first X valid logical index
 */
#define STARTINGX(v) ((v).xinit)

/** Returns the last X valid logical index
 */
#define FINISHINGX(v) ((v).xinit + (int)(v).xdim - 1)

/** Returns the first Y valid logical index
 */
#define STARTINGY(v) ((v).yinit)

/** Returns the last Y valid logical index
 */
#define FINISHINGY(v) ((v).yinit + (int)(v).ydim - 1)

/** Returns the first Z valid logical index
 */
#define STARTINGZ(v) ((v).zinit)

/** Returns the last Z valid logical index
 */
#define FINISHINGZ(v) ((v).zinit + (int)(v).zdim - 1)

/** Check if x is inside logical bounds
 */
#define INSIDEX(v, x) ((x) >= STARTINGX(v) && (x) <= FINISHINGX(v))
/** Check if y is inside logical bounds
 */
#define INSIDEY(v, y) ((y) >= STARTINGY(v) && (y) <= FINISHINGY(v))
/** Check if z is inside logical bounds
 */
#define INSIDEZ(v, z) ((z) >= STARTINGZ(v) && (z) <= FINISHINGZ(v))
/** Check if a position x, y is inside the logical index bounds
 */
#define INSIDEXY(v, x, y) (INSIDEX(v, x) && INSIDEY(v, y))
/** Check if a position x, y is inside the logical index bounds
 */
#define INSIDEXYZ(v, x, y, z) (INSIDEX(v, x) && INSIDEY(v, y) && INSIDEZ(v,z))

/** Access to X dimension (size)
 */
#define XSIZE(v) ((v).xdim)

/** Access to Y dimension (size)
 */
#define YSIZE(v) ((v).ydim)

/** Access to Z dimension (size)
 */
#define ZSIZE(v) ((v).zdim)

/** Access to N dimension (size)
 */
#define NSIZE(v) ((v).ndim)

/** Access to XY dimension (Ysize*Xsize)
 */
#define YXSIZE(v) ((v).yxdim)

/** Access to XYZ dimension (Zsize*Ysize*Xsize)
 */
#define ZYXSIZE(v) ((v).zyxdim)

/** Access to XYZN dimension (Nsize*Zsize*Ysize*Xsize)
 */
#define MULTIDIM_SIZE(v) ((v).nzyxdim)

/** Access to XYZN dimension (Nsize*Zsize*Ysize*Xsize)
 */
#define NZYXSIZE(v) ((v).nzyxdim)

/** Array access.
 *
 * This macro gives you access to the array (T **)
 */
#ifndef MULTIDIM_ARRAY
#define MULTIDIM_ARRAY(v) ((v).data)
#endif

/** Access to a direct element.
 * v is the array, l is the image, k is the slice, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_NZYX_ELEM(v, l, k, i, j) ((v).data[(l)*ZYXSIZE(v)+(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Access to a direct element.
 * v is the array, l is the image, k is the slice, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_ZYX_ELEM(v, k, i, j) ((v).data[(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Access to a direct element.
 * v is the array, l is the image, k =0, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_N_YX_ELEM(v, l, i, j) ((v).data[(l)*ZYXSIZE(v)             +((i)*XSIZE(v))+(j)])
/** Access to a direct element.
 * v is the array, l is the image, k =0, i = 0 and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_N__X_ELEM(v, l, j) ((v).data[(l)*ZYXSIZE(v)+(j)])

/** Multidim element: Logical access.
 */
#define NZYX_ELEM(v, l, k, i, j)  \
    DIRECT_NZYX_ELEM((v), (l), (k) - STARTINGZ(v), (i) - STARTINGY(v), (j) - STARTINGX(v))

/** Access to a direct element.
 * v is the array, k is the slice and n is the number of the pixel (combined i and j)
 * within the slice.
 */
#define DIRECT_MULTIDIM_ELEM(v,n) ((v).data[(n)])

/** For all direct elements in the array
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'n' which goes over the slices and 'n' that
 * goes over the pixels in each slice.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
 * {
 *     std::cout << DIRECT_MULTIDIM_ELEM(v,n) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v) \
    for (size_t n=0; n<NZYXSIZE(v); ++n)

/** For all direct elements in the array
 *
 * This macro is used to generate loops for the array in an easy
 * manner. It defines internal indexes 'l', 'k','i' and 'j' which
 * ranges over the n volume using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v)
 * {
 *     std::cout << DIRECT_NZYX_ELEM(v,l, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V) \
    for (size_t l=0; l<NSIZE(V); ++l) \
        for (size_t k=0; k<ZSIZE(V); ++k) \
            for (size_t i=0; i<YSIZE(V); ++i)      \
                for (size_t j=0; j<XSIZE(V); ++j)

/** For all direct elements in the array
 *
 * This macro is used to generate loops for the array in an easy
 * manner. It defines internal indexes 'l', 'k','i' and 'j' which
 * ranges over the n volume using its logical definition.
 *
 * @code
 * FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v)
 * {
 *     std::cout << NZYX_ELEM(v,l, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V) \
    for (size_t l=0; l<NSIZE(V); ++l) \
        for (int k=STARTINGZ(V); k<=FINISHINGZ(V); ++k) \
            for (int i=STARTINGY(V); i<=FINISHINGY(V); ++i)     \
                for (int j=STARTINGX(V); j<=FINISHINGX(V); ++j)

/** For all direct elements in the array, pointer version
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'k' which goes over the slices and 'n' that
 * goes over the pixels in each slice. Each element can be accessed through
 * an external pointer called ptr.
 *
 * @code
 * T* ptr=NULL;
 * size_t n;
 * FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
 * {
 *     std::cout << *ptr << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr) \
    for ((n)=0, (ptr)=(v).data; (n)<NZYXSIZE(v); ++(n), ++(ptr))

/** Access to a direct element.
 * v is the array, k is the slice (Z), i is the Y index and j is the X index.
 */
#define DIRECT_A3D_ELEM(v,k,i,j) ((v).data[YXSIZE(v)*(k)+((i)*XSIZE(v))+(j)])

/** A short alias for the previous function.
 *
 */
#define dAkij(V, k, i, j) DIRECT_A3D_ELEM(V, k, i, j)

/** Volume element: Logical access.
 *
 * @code
 * A3D_ELEM(V, -1, -2, 1) = 1;
 * val = A3D_ELEM(V, -1, -2, 1);
 * @endcode
 */
#define A3D_ELEM(V, k, i, j) \
    DIRECT_A3D_ELEM((V),(k) - STARTINGZ(V), (i) - STARTINGY(V), (j) - STARTINGX(V))

/** For all elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(V)
 * {
 *     std::cout << V(k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY3D(V) \
    for (int k=STARTINGZ(V); k<=FINISHINGZ(V); ++k) \
        for (int i=STARTINGY(V); i<=FINISHINGY(V); ++i) \
            for (int j=STARTINGX(V); j<=FINISHINGX(V); ++j)

/** For all elements in common.
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two volumes in an easy manner. Then k, i and j (locally defined)
 * range from
 *
 * MAX(STARTINGZ(V1),STARTINGZ(V2)) to MIN(FINISHINGZ(V1),FINISHINGZ(V2)),
 * MAX(STARTINGY(V1),STARTINGY(V2)) to MIN(FINISHINGY(V1),FINISHINGY(V2)),
 * MAX(STARTINGX(V1),STARTINGX(V2)) to MIN(FINISHINGX(V1),FINISHINGX(V2))
 *
 * (included limits) respectively. You need to define SPEED_UP_temps.
 *
 * @code
 * SPEED_UP_temps;
 * MultidimArray< double > V1(10, 10, 10), V2(20, 20, 20);
 * V1.setXmippOrigin();
 * V2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(V1, V2)
 * {
 *    // ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(V1, V2) \
    ispduptmp0 = XMIPP_MAX(STARTINGZ(V1), STARTINGZ(V2)); \
    ispduptmp1 = XMIPP_MIN(FINISHINGZ(V1),FINISHINGZ(V2)); \
    ispduptmp2 = XMIPP_MAX(STARTINGY(V1), STARTINGY(V2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(V1),FINISHINGY(V2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(V1), STARTINGX(V2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(V1),FINISHINGX(V2)); \
    for (int k=ispduptmp0; k<=ispduptmp1; ++k) \
        for (int i=ispduptmp2; i<=ispduptmp3; ++i) \
            for (int j=ispduptmp4; j<=ispduptmp5; ++j)

/** For all direct elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V)
 * {
 *     std::cout << DIRECT_A3D_ELEM(m, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V) \
    for (size_t k=0; k<ZSIZE(V); ++k) \
        for (size_t i=0; i<YSIZE(V); ++i) \
            for (size_t j=0; j<XSIZE(V); ++j)

/** Access to a direct element of a matrix.
 * v is the array, i and j define the element v_ij.
 *
 * Be careful because this is physical access, usually matrices follow the C
 * convention of starting index==0 (X and Y). This function should not be used
 * as it goes against the vector library philosophy unless you explicitly want
 * to access directly to any value in the matrix without taking into account its
 * logical position
 *
 * @code
 * DIRECT_A2D_ELEM(m, 0, 0) = 1;
 * val = DIRECT_A2D_ELEM(m, 0, 0);
 * @endcode
 */
#define DIRECT_A2D_ELEM(v,i,j) ((v).data[(i)*(v).xdim+(j)])

/** Short alias for DIRECT_A2D_ELEM
 */
#define dAij(M, i, j) DIRECT_A2D_ELEM(M, i, j)

/** Matrix element: Logical access
 *
 * @code
 * A2D_ELEM(m, -2, 1) = 1;
 * val = A2D_ELEM(m, -2, 1);
 * @endcode
 */
#define A2D_ELEM(v, i, j) \
    DIRECT_A2D_ELEM(v, (i) - STARTINGY(v), (j) - STARTINGX(v))

/** TRUE if both arrays have the same shape
 *
 * Two arrays have the same shape if they have the same size and the same
 * starting point. Be aware that this is a macro which simplifies to a boolean.
 */
#define SAME_SHAPE2D(v1, v2) \
    (XSIZE(v1) == XSIZE(v2) && \
     YSIZE(v1) == YSIZE(v2) && \
     STARTINGX(v1) == STARTINGX(v2) && \
     STARTINGY(v1) == STARTINGY(v2))

#define SAME_SHAPE3D(v1, v2) \
    (XSIZE(v1) == XSIZE(v2) && \
     YSIZE(v1) == YSIZE(v2) && \
     ZSIZE(v1) == ZSIZE(v2) && \
     STARTINGX(v1) == STARTINGX(v2) && \
     STARTINGY(v1) == STARTINGY(v2) && \
     STARTINGZ(v1) == STARTINGZ(v2))


/** For all elements in the array
 *
 * This macro is used to generate loops for the matrix in an easy way. It
 * defines internal indexes 'i' and 'j' which ranges the matrix using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY2D(m)
 * {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY2D(m) \
    for (int i=STARTINGY(m); i<=FINISHINGY(m); ++i) \
        for (int j=STARTINGX(m); j<=FINISHINGX(m); ++j)

/** For all elements in common
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two images in an easy manner. Then i and j (locally defined) range
 * from MAX(STARTINGY(V1), STARTINGY(V2)) to MIN(FINISHINGY(V1),
 * FINISHINGY(V2)), MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1),
 * FINISHINGX(V2)) (included limits) respectively. You need to define
 * SPEED_UP_temps.
 *
 * @code
 * MultidimArray< double > m1(10, 10), m2(20, 20);
 * m1.setXmippOrigin();
 * m2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY2D(m1, m2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY2D(m1, m2) \
    ispduptmp2 = XMIPP_MAX(STARTINGY(m1), STARTINGY(m2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(m1), FINISHINGY(m2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(m1), STARTINGX(m2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(m1), FINISHINGX(m2)); \
    for (int i=ispduptmp2; i<=ispduptmp3; ++i) \
        for (int j=ispduptmp4; j<=ispduptmp5; ++j)

/** For all elements in the array between corners.
 *
 * This macro is used to generate loops for a volume in an easy manner. Then
 *  ZZ(r), YY(r) and XX(r) range from
 *
 * (int) ZZ(corner1) to (int)ZZ(corner2),
 * (int) YY(corner1) to (int)YY(corner2),
 * (int) XX(corner1) to (int) XX(corner2) (included limits) respectively.
 *
 * Notice that corner1 and corner2 need only be MultidimArray.
 *
 * @code
 * MultidimArray< double > corner1(3), corner2(3), r(3);
 * XX(corner1) = -1; XX(corner2) = 1;
 * YY(corner1) = -2; YY(corner2) = 2;
 * ZZ(corner1) = -3; ZZ(corner2) = 3;
 *
 * FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(corner1, corner2) \
    for (ZZ(r)=ZZ((corner1)); ZZ(r)<=ZZ((corner2)); ++ZZ(r)) \
        for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); ++YY(r)) \
            for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); ++XX(r))

/** For all elements in the array between corners
 *
 * This macro is used to generate loops for a matrix in an easy manner. It needs
 * an externally defined MultidimArray< double > r(2). Then YY(r) and XX(r) range
 * from (int) YY(corner1) to (int)YY(corner2), (int) XX(corner1) to (int)
 * XX(corner2) (included limits) respectively. Notice that corner1 and corner2
 * need only be MultidimArray.
 *
 * @code
 * MultidimArray< double > corner1(2), corner2(2);
 * MultidimArray< int > r(2);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * YY(corner1) = -2;
 * YY(corner2) = 2;
 *
 * FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2) \
    for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); ++YY(r)) \
        for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); ++XX(r))

/** For all elements in the array between corners
 *
 * This macro is used to generate loops for a vector in an easy manner. It needs
 * an externally defined MultidimArray< double > r(1). Then XX(r) ranges from
 * (int) XX(corner1) to (int) XX(corner2) (included limits) (notice that corner1
 * and corner2 need only to be MultidimArray).
 *
 * @code
 * MultidimArray< double > corner1(1), corner2(1), r(1);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * FOR_ALL_ELEMENTS_IN_ARRAY1D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(XX(r)) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY1D_BETWEEN(corner1, corner2) \
    for (XX(r)=(int) XX((corner1)); XX(r)<=(int) XX((corner2)); ++XX(r))

/** For all elements in the array, accessed physically
 *
 * This macro is used to generate loops for the matrix in an easy way using
 * physical indexes. It defines internal indexes 'i' and 'j' which ranges the
 * matrix using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(m)
 * {
 *     std::cout << DIRECT_A2D_ELEM(m, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(m) \
    for (size_t i=0; i<YSIZE(m); ++i) \
        for (size_t j=0; j<XSIZE(m); ++j)

/** Vector element: Physical access
 *
 * Be careful because this is physical access, usually vectors follow the C
 * convention of starting index==0. This function should not be used as it goes
 * against the vector library philosophy unless you explicitly want to access
 * directly to any value in the vector without taking into account its logical
 * position.
 *
 * @code
 * DIRECT_A1D_ELEM(v, 0) = 1;
 * val = DIRECT_A1D_ELEM(v, 0);
 * @endcode
 */
#define DIRECT_A1D_ELEM(v, i) ((v).data[(i)])

/** A short alias to previous function
 */
#define dAi(v, i) DIRECT_A1D_ELEM(v, i)

/** Vector element: Logical access
 *
 * @code
 * A1D_ELEM(v, -2) = 1;
 * val = A1D_ELEM(v, -2);
 * @endcode
 */
#define A1D_ELEM(v, i) DIRECT_A1D_ELEM(v, (i) - ((v).xinit))

/** For all elements in the array
 *
 * This macro is used to generate loops for the vector in an easy manner. It
 * defines an internal index 'i' which ranges the vector using its mathematical
 * definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY1D(v)
 * {
 *     std::cout << v(i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY1D(v) \
    for (int i=STARTINGX(v); i<=FINISHINGX(v); ++i)

/** For all elements in common
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two vectors in an easy manner. Then i (locally defined) ranges from
 * MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1), FINISHINGX(V2))
 * (included limits) respectively. You need to define SPEED_UP_temps.
 *
 * @code
 * MultidimArray< double > v1(10), v2(20);
 * v1.setXmippOrigin();
 * v2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY1D(v1, v2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY1D(v1, v2) \
    ispduptmp4 = XMIPP_MAX(STARTINGX(v1), STARTINGX(v2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(v1), FINISHINGX(v2)); \
    for (int i=ispduptmp4; i<=ispduptmp5; ++i)

/** For all elements in the array, accessed physically
 *
 * This macro is used to generate loops for the vector in an easy way using
 * physical indexes. It defines internal the index 'i' which ranges the vector
 * using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v)
 * {
 *     std::cout << DIRECT_A2D_ELEM(v, i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v) \
    for (size_t i=0; i<v.xdim; ++i)

/** Macro to check whether a point is inside or outside a given matrix. */
#define OUTSIDE(i,j) \
            ((j) < STARTINGX(*this) || (j) > FINISHINGX(*this) || \
             (i) < STARTINGY(*this) || (i) > FINISHINGY(*this))

/** Macro to check whether a point is inside or outside a given matrix. */
#define OUTSIDE3D(k, i,j) \
            ((j) < STARTINGX(*this) || (j) > FINISHINGX(*this) || \
             (i) < STARTINGY(*this) || (i) > FINISHINGY(*this) || \
             (k) < STARTINGZ(*this) || (k) > FINISHINGZ(*this))
//@}

// Look up table lenght to be used in interpolation.
#define     LOOKUP_TABLE_LEN        6


/** Template class for Xmipp arrays.
  * This class provides physical and logical access.
*/
class MultidimArrayBase
{
public:
    // Destroy data
    bool destroyData;

    // Number of images
    size_t ndim;

    // Number of elements in Z
    size_t zdim;

    // Number of elements in Y
    size_t ydim;

    // Number of elements in X
    size_t xdim;

    // Number of elements in YX
    size_t yxdim;

    // Number of elements in ZYX
    size_t zyxdim;

    // Number of elements in NZYX
    size_t nzyxdim;

    // Z init
    int zinit;

    // Y init
    int yinit;

    // X init
    int xinit;

    //Alloc memory or map to a file
    bool     mmapOn;
    //Mapped file handler
    FILE*      mFd;
    // Number of elements in NZYX in allocated memory
    size_t nzyxdimAlloc;
public:
    virtual ~MultidimArrayBase()
    {}

    // Virtual declarations to be used from MultidimArrayGeneric
    virtual void clear() = 0;
    virtual void selfReverseX() = 0;
    virtual void selfReverseY() = 0;
    virtual void selfReverseZ() = 0;
    virtual double computeAvg() const = 0;
    virtual void computeDoubleMinMaxRange(double& minval, double& maxval,size_t offset, size_t size) const = 0;
    virtual void maxIndex(size_t &lmax, int& kmax, int& imax, int& jmax) const = 0;
    virtual void coreAllocateReuse() = 0;
    virtual void coreDeallocate()= 0;

    /* return the value of the data pointer
     */
    virtual void * getArrayPointer() const = 0;

    /// @name Size
    //@{

    /** Sets new N dimension.
        *
        *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
        */
    void setNdim(int Ndim);

    /** Sets new Z dimension.
     *
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     *
     */
    void setZdim(int Zdim);

    /** Sets new Y dimension.
     *
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     *
     */
    void setYdim(int Ydim);

    /** Sets new X dimension.
      *
      *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
      *
      */
    void setXdim(int Xdim);

    /** Sets new 4D dimensions.
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     */
    void setDimensions(int Xdim, int Ydim, int Zdim, size_t Ndim);

    /** Sets new 4D dimensions.
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     */
    void setDimensions(ArrayDim &newDim);

    /** Get the array dimensions.
     */
    void getDimensions(size_t& Xdim, size_t& Ydim, size_t& Zdim, size_t &Ndim) const;
    void getDimensions(ArrayDim &idim) const;
    ArrayDim getDimensions() const;


    /** Get dimensions.
     *
     * Returns the size of the object in a 4D vector. If the object is a matrix
     * or a vector, then the higher order dimensions will be set to 1, ie,
     * (Xdim, 1, 1) or (Xdim, Ydim, 1).
     *
     * This function is not ported to Python.
     */
    void getDimensions(int* size) const;

    /** Returns the total size of the multidimArray
     *
     * @code
     * if (V.getSize() > 1) ...
     * @endcode
     */
    size_t getSize() const;

    /** Resize to a given size
      *
      * This function resize the actual array to the given size. The origin is
      * not modified. If the actual array is larger than the pattern then the
      * values outside the new size are lost, if it is smaller then 0's are
      * added. An exception is thrown if there is no memory.
      *
      * @code
      * V1.resize(3, 3, 2);
      * @endcode
      */
    virtual void resize(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim, bool copy=true) = 0;

    /** Resize a single 3D image
     *
     * This function assumes n is 1
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
     */
    void resize(size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(1, Zdim, Ydim, Xdim);
    }

    /** Resize a single 2D image
     *
     * This function assumes n and z are 1
     * @code
     * V1.resize(3, 2);
     * @endcode
     */
    void resize(size_t Ydim, size_t Xdim)
    {
        resize(1, 1, Ydim, Xdim);
    }

    /** Resize a single 1D image
     *
     * This function assumes n and z and y are 1
     * @code
     * V1.resize(2);
     * @endcode
     */
    void resize(size_t Xdim)
    {
        resize(1, 1, 1, Xdim);
    }

    /** Resize an image using the dimensions
     *  from an ArrayDim structure.
     */
    void resize(ArrayDim &adim, bool copy=true);

    /** Resize with no copy a single 3D image
        */
    void resizeNoCopy(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(Ndim, Zdim, Ydim, Xdim, false);
    }

    /** Resize with no copy a single 3D image
        */
    void resizeNoCopy(size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(1, Zdim, Ydim, Xdim, false);
    }

    /** Resize a single 2D image with no copy
     */
    void resizeNoCopy(size_t Ydim, size_t Xdim)
    {
        resize(1, 1, Ydim, Xdim, false);
    }

    /** Resize a single 1D image with no copy
     */
    void resizeNoCopy(size_t Xdim)
    {
        resize(1, 1, 1, Xdim, false);
    }

    /** Returns Y dimension.
       */
    inline size_t rowNumber() const
    {
        return ydim;
    }

    /** Returns X dimension.
     */
    inline size_t colNumber() const
    {
        return xdim;
    }

    /** Copy the shape parameters
      *
      */
    void copyShape(const MultidimArrayBase &m);

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    inline bool sameShape(const MultidimArrayBase &op) const
    {
        return (NSIZE(*this) == NSIZE(op) &&
                XSIZE(*this) == XSIZE(op) &&
                YSIZE(*this) == YSIZE(op) &&
                ZSIZE(*this) == ZSIZE(op) &&
                STARTINGX(*this) == STARTINGX(op) &&
                STARTINGY(*this) == STARTINGY(op) &&
                STARTINGZ(*this) == STARTINGZ(op));
    }

    /** Set logical origin in Xmipp fashion.
      *
      * This function adjust the starting points in the array such that the
      * center of the array is defined in the Xmipp fashion.
      *
      * @code
      * V.setXmippOrigin();
      * @endcode
      */
    void setXmippOrigin();

    /** Reset logical origin to zeros.
     *
     * This function adjust the starting points in the array such
     * that upper left corner begins in zero.
     *
     * @code
     * V.resetOrigin();
     * @endcode
     */
    void resetOrigin();

    /** Move origin to.
      *
      * This function adjust logical indexes such that the Xmipp origin of the
      * array moves to the specified position. For instance, an array whose x
      * indexes go from -1 to 1, if we move the origin to 4, then the x indexes
      * go from 3 to 5. This is very useful for convolution operations where you
      * only need to move the logical starting of the array.
      *
      */
    void moveOriginTo(int k, int i, int j);

    /** Move origin to.
      *
      * This function adjust logical indexes such that the Xmipp origin of the
      * array moves to the specified position. For instance, an array whose x
      * indexes go from -1 to 1, if we move the origin to 4, then the x indexes
      * go from 3 to 5. This is very useful for convolution operations where you
      * only need to move the logical starting of the array.
      *
      */
    void moveOriginTo(int i, int j);

    /** Returns the first valid logical Z index.
      */
    inline int startingZ() const
    {
        return zinit;
    }

    /** Returns the last valid logical Z index.
     */
    inline int finishingZ() const
    {
        return zinit + zdim - 1;
    }

    /** Returns the first valid logical Y index.
     */
    inline int startingY() const
    {
        return yinit;
    }

    /** Returns the last valid logical Y index.
     */
    inline int finishingY() const
    {
        return yinit + ydim - 1;
    }

    /** Returns the first valid logical X index.
     */
    inline int startingX() const
    {
        return xinit;
    }

    /** Returns the last valid logical X index.
     */
    inline int finishingX() const
    {
        return xinit + xdim - 1;
    }

    /** IsCorner (in 2D or 3D matrix)
         *
         * TRUE if the logical index given is a corner of the definition region of this
         * array.
         */
    bool isCorner(const Matrix1D< double >& v) const;

    /** Outside for 3D matrices
       *
       * TRUE if the logical index given is outside the definition region of this
       * array.
       */
    inline bool outside(int k, int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this) ||
                k < STARTINGZ(*this) || k > FINISHINGZ(*this));
    }

    /** Outside for 2D matrices
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    inline bool outside(int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this));
    }

    /** Outside for 1D matrices
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    inline bool outside(int i) const
    {
        return (i < STARTINGX(*this) || i > FINISHINGX(*this));
    }

    /** Outside
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(const Matrix1D<double> &r) const;

    //@}


    /** Returns the multidimArray dimension.
     *
     * @code
     * int dim = V.getDim();
     * @endcode
     */
    inline int getDim() const
    {
        if (NZYXSIZE(*this) < 1)
            return 0;
        if (NSIZE(*this) > 1)
            return 4;
        if (ZSIZE(*this) > 1)
            return 3;
        if (YSIZE(*this) > 1)
            return 2;
        return 1;
    }

    /** Sets mmap.
     *
     * Sets on/off mmap flag to allocate memory in a file.
     *
     */
    void setMmap(bool mmap)
    {
        coreDeallocate();
        mmapOn = mmap;
    }

    void maxIndex(ArrayCoord &pos) const
    {
        maxIndex(pos.n, pos.z, pos.y, pos.x);
    }

    /** 3D Indices for the maximum element.
      *
      * This function just calls to the 4D function
      */
    void maxIndex(int& kmax, int& imax, int& jmax) const
    {
        size_t dum;
        maxIndex(dum, kmax, imax, jmax);
    }

    /** 2D Indices for the maximum element.
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& imax, int& jmax) const
    {
        size_t dum;
        int idum;
        maxIndex(dum, idum, imax, jmax);
    }

    /** 1D Indices for the maximum element.
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& jmax) const
    {
        size_t dum;
        int idum;
        maxIndex(dum, idum, idum, jmax);
    }

    /** Print shape of multidimensional array.
     *
     * This function shows the size, starting and finishing indexes of the
     * given array. No end of line is printed neither at the beginning nor
     * the end.
     *
     * @code
     * v.printShape();
     *
     * std::ofstream fh;
     * ...;
     * v.printShape(fh);
     * @endcode
     */
    void printShape(std::ostream& out = std::cout) const;
};

#endif /* XMIPPCORE_CORE_MULTIDIM_ARRAY_BASE_H_ */
