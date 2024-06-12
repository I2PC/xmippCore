/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
 *    Joaquin Oton       (joton@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#ifndef CORE_IMAGE_H
#define CORE_IMAGE_H

#include <typeinfo>
#include <set>
#include "multidim_array.h"
#include "xmipp_image_base.h"
#include "xmipp_memory.h"
#include "utils/half.hpp"

/// @addtogroup Images
//@{

/** Size of the page used to read and write images from/to file */
const size_t rw_max_page_size = 4194304; // 4Mb

/** Template class for images.
 * The image class is the general image handling class.
 */
template<typename T>
class Image : public ImageBase
{

public:
    MultidimArray<T> data; // The image data array

public:
    /** Empty constructor
     *
     * An empty image is created.
     *
     * @code
     * Image<double> I;
     * @endcode
     */
    Image()
    {
        mdaBase = (MultidimArrayBase*) &data;
        init();
    }

    /** Constructor with size and filename
     *
     * An image file, which name and format are given by filename,
     * is created with the given size. Then the image is mapped to this file.
     *
     * @code
     * Image I(64,64,1,1,"image.spi");
     * @endcode
     */
    Image(int Xdim, int Ydim, int Zdim, int Ndim, const FileName &_filename)
    {
        mdaBase = (MultidimArrayBase*) &data;
        init();
        mmapOnWrite = true;
        data.setDimensions(Xdim, Ydim, Zdim, Ndim);
        MD.resize(Ndim);
        filename = _filename;
        ImageFHandler *hFile = openFile(_filename, WRITE_OVERWRITE);
        _write(_filename, hFile, ALL_IMAGES, false, WRITE_OVERWRITE);
        closeFile(hFile);
    }

    /** Constructor with size
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X. If _mmapOn is True then image is allocated
     * in a temporary file.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    Image(int Xdim, int Ydim, int Zdim = 1, int Ndim = 1, bool _mmapOn = false)
    {
        mdaBase = (MultidimArrayBase*) &data;
        init();
        data.setMmap(_mmapOn);
        data.coreAllocate(Ndim, Zdim, Ydim, Xdim);
        MD.resize(Ndim);
    }

    Image(const Image<T> &im)
    {
        mdaBase = (MultidimArrayBase*) &data;
        init();
        *this = im;
    }

    /** Constructor with MultidimArray alias
     *
     *  An image is created directly with its multidimarray aliased to im.
     *  This function is useful when debugging and you want to save multidimarrays.
     */
    Image(const MultidimArray<T> &im)
    {
        mdaBase = (MultidimArrayBase*) &data;
        init();
        data.alias(im);
    }

    /** Destructor.
     */
    virtual
    ~Image()
    {
        clearData();
    }

    /** Clear.
     * Initialize everything to 0
     */
    void
    clear()
    {
        clearData();
        init();
    }

    /** Clear the data pointer either mapped or not.
     *  Initially fix the pointer position in case of
     *  tricking a lower number of slices and starting
     *  in a slice different from zero.
     */
    void
    clearData()
    {
        // If MultidimArray pointer has been moved to a slice different from zero, then reset it.
        // This check must be done prior to mappedSize check, since mappedSlice is a trick over data pointer
        if (virtualOffset != 0)
            movePointerTo(ALL_SLICES);
        if (mmapOnRead || mmapOnWrite)
            munmapFile();

        else
            data.clear();

    }

    /** Check whether image is complex based on T
     */
    bool
    isComplexT() const
    {
        return (typeid(T) == typeid(std::complex<double>)
                || typeid(T) == typeid(std::complex<float>));
    }


    /**
     * Obtain the inverse axis mapping
    */
    void
    getInverseAxisOrder(const std::array<int,4> &order, 
                        std::array<int,4> &result )
    {
        for (size_t i = 0; i < result.size(); ++i)
        {
            size_t j = 0;
            while(j < order.size() && order[j] != i) ++j; // Find inverse mapping
            if (j >= order.size())
                REPORT_ERROR(ERR_LOGIC_ERROR, "Invalid axis mapping");
            
            result[i] = j;
        }
    }

    /**
     * Trasposes the given size array according to inverseOrder
    */
    void
    transposeAxisSizes(const std::array<size_t,4> &sizes, 
                       const std::array<int,4> &order,
                       std::array<size_t,4> &result )
    {
        std::array<int, 4> inverseOrder;
        getInverseAxisOrder(order, inverseOrder);

        result = {
            sizes[inverseOrder[0]],
            sizes[inverseOrder[1]],
            sizes[inverseOrder[2]],
            sizes[inverseOrder[3]]
        };
    }

    /**
     * Trasposes the given MultidimArray with the given order
    */
   void
   transposeInPlace(MultidimArray<T> &multidimArray, const std::array<int,4> &order)
   {
        std::array<int, 4> inverseOrder;
        getInverseAxisOrder(order, inverseOrder);

        // Creating new multidim array of the same size than the original
        const std::array<size_t,4> sizes = {
            NSIZE(multidimArray),
            ZSIZE(multidimArray),
            YSIZE(multidimArray),
            XSIZE(multidimArray)
        };

        MultidimArray<T> result(
            sizes[inverseOrder[0]],
            sizes[inverseOrder[1]],
            sizes[inverseOrder[2]],
            sizes[inverseOrder[3]]
        );

        // Performing transposition in a loop for every dimension
        for (size_t n = 0; n < NSIZE(multidimArray); n++) {
            for (size_t z = 0; z < ZSIZE(multidimArray); z++) {
                for (size_t y = 0; y < YSIZE(multidimArray); y++) {
                    for (size_t x = 0; x < XSIZE(multidimArray); x++) {
                        // Defining array to access with the axis orders
                        const std::array<size_t,4> indices = {n, z, y, x};
                        const auto l = indices[inverseOrder[0]];
                        const auto k = indices[inverseOrder[1]];
                        const auto i = indices[inverseOrder[2]];
                        const auto j = indices[inverseOrder[3]];

                        // Transposing element
                        DIRECT_NZYX_ELEM(result, l, k, i, j) = DIRECT_NZYX_ELEM(multidimArray, n, z, y, x);
                    }
                }
            }
        }

        // Remapping pointers from original multidim array to transposed one
        multidimArray = std::move(result);
   }

    /** Cast a page of data from type dataType to type Tdest
     *    input pointer  char *
     */
    void
    castPage2T(char * page, T * ptrDest, DataType datatype, size_t pageSize)
    {
        switch (datatype)
        {
            case DT_Unknown:
                REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR: datatype is Unknown_Type");
            case DT_UHalfByte:
            case DT_UChar:
            {
                if (typeid(T) == typeid(unsigned char))
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                else
                {
                    const auto* ptr = (unsigned char *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_SChar:
            {
                if (typeid(T) == typeid(signed char))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (signed char *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_UShort:
            {
                if (typeid(T) == typeid(unsigned short))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (unsigned short *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_Short:
            {
                if (typeid(T) == typeid(short))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (short *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_UInt:
            {
                if (typeid(T) == typeid(unsigned int))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (unsigned int *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_Int:
            {
                if (typeid(T) == typeid(int))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (int *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_Long:
            {
                if (typeid(T) == typeid(long))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (long *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_Float:
            {
                if (typeid(T) == typeid(float))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (float *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_Double:
            {
                if (typeid(T) == typeid(double))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (double *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            case DT_HalfFloat:
            {
                if (typeid(T) == typeid(half_float::half))
                {
                    memcpy(ptrDest, page, pageSize * sizeof(T));
                }
                else
                {
                    const auto* ptr = (half_float::half *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        ptrDest[i] = (T) *ptr;
                }
                break;
            }
            default:
            {
                std::cerr << "Datatype= " << datatype << std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT, " ERROR: cannot cast datatype to T");
            }
        }

    }

    /** Cast page from T to datatype
     *  input pointer char *
     */
    void
    castPage2Datatype(T * srcPtr, char * page, DataType datatype,
                      size_t pageSize) const
    {
        switch (datatype)
        {
            case DT_Float:
            {
                if (typeid(T) == typeid(float))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (float *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (float) srcPtr[i];
                }
                break;
            }
            case DT_Double:
            {
                if (typeid(T) == typeid(double))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (double *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (double) srcPtr[i];
                }
                break;
            }
            case DT_UShort:
            {
                if (typeid(T) == typeid(unsigned short))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (unsigned short *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (unsigned short) srcPtr[i];
                }
                break;
            }
            case DT_Short:
            {
                if (typeid(T) == typeid(short))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (short *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (short) srcPtr[i];
                }
                break;
            }
            case DT_UHalfByte:
            case DT_UChar:
            {
                if (typeid(T) == typeid(unsigned char))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (unsigned char *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (unsigned char) srcPtr[i];
                }
                break;
            }
            case DT_SChar:
            {
                if (typeid(T) == typeid(char))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (char *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (char) srcPtr[i];
                }
                break;
            }
            case DT_HalfFloat:
            {
                if (typeid(T) == typeid(half_float::half))
                {
                    memcpy(page, srcPtr, pageSize * sizeof(T));
                }
                else
                {
                    auto* ptr = (half_float::half *) page;
                    for (size_t i = 0; i < pageSize; ++i, ++ptr)
                        *ptr = (half_float::half) srcPtr[i];
                }
                break;
            }
            default:
            {
                std::cerr << "outputDatatype = " << datatype << std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT,
                             " ERROR: cannot cast T to outputDatatype");
            }
        }
    }

    /* Convert the pixels values from one datatype to another, taking into account for datatypes
     * of same bitdepth the shift of the minimum values. In other cases, the conversion is done
     * adjusting the input values in the range of output datatype.
     */
    void
    castConvertPage2Datatype(T * srcPtr, char * page, DataType datatype,
                             size_t pageSize, double min0, double max0, CastWriteMode castMode =
    CW_CONVERT) const
    {

        double minF, maxF;
        double slope;
        size_t n;
        DataType myTypeId = myT();

        switch (datatype)
        {
            case DT_UHalfByte:
            case DT_UChar:
            {
                if (castMode == CW_CONVERT && myTypeId == DT_SChar)
                {
                    slope = 1;
                    min0 -= CHAR_MIN;
                }
                else if (castMode == CW_CAST)
                {
                	minF  = 0;
                	min0  = 0;
                	slope = 1;
                }
                else
                {
                    minF = 0;
                    maxF = UCHAR_MAX;
                    if (max0 != min0)
                        slope = static_cast<double>(maxF - minF)
                                / static_cast<double>(max0 - min0);
                    else
                        slope = 0;
                }
                unsigned char * ptr = (unsigned char *) page;

                for (n = 0; n < pageSize; n++)
                    ptr[n] = static_cast<unsigned char>(minF
                                                        + (slope * static_cast<double>(srcPtr[n] - min0)));

                break;
            }
            case DT_SChar:
            {
                if (castMode == CW_CONVERT && myTypeId == DT_UChar)
                {
                    slope = 1;
                    min0 += CHAR_MIN;
                }
                else if (castMode == CW_CAST)
                {
                	minF  = 0;
                	min0  = 0;
                	slope = 1;
                }
                else
                {
                    minF = CHAR_MIN;
                    maxF = CHAR_MAX;
                    if (max0 != min0)
                        slope = static_cast<double>(maxF - minF)
                                / static_cast<double>(max0 - min0);
                    else
                        slope = 0;
                }
                char * ptr = (char *) page;

                for (n = 0; n < pageSize; n++)
                    ptr[n] = static_cast<char>(minF
                                               + (slope * static_cast<double>(srcPtr[n] - min0)));

                break;
            }
            case DT_UShort:
            {
                if (castMode == CW_CONVERT
                    && (myTypeId == DT_SChar || myTypeId == DT_Short))
                {
                    slope = 1;
                    min0 -= SHRT_MIN;
                }
                else if (castMode == CW_CONVERT && (myTypeId == DT_UChar))
                {
                    slope = 1;
                }
                else if (castMode == CW_CAST)
                {
                	minF  = 0;
                	min0  = 0;
                	slope = 1;
                }
                else
                {
                    minF = 0;
                    maxF = USHRT_MAX;
                    if (max0 != min0)
                        slope = static_cast<double>(maxF - minF)
                                / static_cast<double>(max0 - min0);
                    else
                        slope = 0;
                }

                unsigned short * ptr = (unsigned short *) page;

                for (n = 0; n < pageSize; n++)
                    ptr[n] = static_cast<unsigned short>(minF
                                                         + (slope * static_cast<double>(srcPtr[n] - min0)));

                break;
            }
            case DT_Short:
            {
                if (castMode == CW_CONVERT
                    && (myTypeId == DT_UChar || myTypeId == DT_UShort))
                {
                    slope = 1;
                    min0 += SHRT_MIN;
                }
                else if (castMode == CW_CONVERT && (myTypeId == DT_SChar))
                {
                    slope = 1;
                }
                else if (castMode == CW_CAST)
                {
                	minF  = 0;
                	min0  = 0;
                	slope = 1;
                }
                else
                {
                    minF = SHRT_MIN;
                    maxF = SHRT_MAX;
                    if (max0 != min0)
                        slope = static_cast<double>(maxF - minF)
                                / static_cast<double>(max0 - min0);
                    else
                        slope = 0;
                }
                short * ptr = (short *) page;

                for (n = 0; n < pageSize; n++)
                    ptr[n] = static_cast<short>(minF
                                                + (slope * static_cast<double>(srcPtr[n] - min0)));

                break;
            }
            case DT_UInt:
            {
                if (castMode == CW_CONVERT
                    && (myTypeId == DT_SChar || myTypeId == DT_Short
                        || myTypeId == DT_Int))
                {
                    slope = 1;
                    min0 -= INT_MIN;
                }
                else if (castMode == CW_CONVERT
                         && (myTypeId == DT_UShort || myTypeId == DT_UChar))
                {
                    slope = 1;
                }
                else if (castMode == CW_CAST)
                {
                	minF  = 0;
                	min0  = 0;
                	slope = 1;
                }
                else
                {
                    minF = 0;
                    maxF = UINT_MAX;
                    if (max0 != min0)
                        slope = static_cast<double>(maxF - minF)
                                / static_cast<double>(max0 - min0);
                    else
                        slope = 0;
                }
                unsigned int * ptr = (unsigned int *) page;

                for (n = 0; n < pageSize; n++)
                    ptr[n] = static_cast<unsigned int>(minF
                                                       + (slope * static_cast<double>(srcPtr[n] - min0)));
                break;
            }
            case DT_Int:
            {
                if (castMode == CW_CONVERT
                    && (myTypeId == DT_UChar || myTypeId == DT_UShort
                        || myTypeId == DT_UInt))
                {
                    slope = 1;
                    min0 += INT_MIN;
                }
                else if (castMode == CW_CONVERT
                         && (myTypeId == DT_Short || myTypeId == DT_SChar))
                {
                    slope = 1;
                }
                else if (castMode == CW_CAST)
                {
                	minF  = 0;
                	min0  = 0;
                	slope = 1;
                }
                else
                {
                    minF = INT_MIN;
                    maxF = INT_MAX;
                    if (max0 != min0)
                        slope = static_cast<double>(maxF - minF)
                                / static_cast<double>(max0 - min0);
                    else
                        slope = 0;
                }
                int * ptr = (int *) page;

                for (n = 0; n < pageSize; n++)
                    ptr[n] = static_cast<int>(minF
                                              + (slope * static_cast<double>(srcPtr[n] - min0)));
                break;
            }
            default:
                castPage2Datatype(srcPtr, page, datatype, pageSize);
        }

    }

    void
    setPage2T(size_t offset, char * page, DataType datatype, size_t pageSize)
    {
        castPage2T(page, MULTIDIM_ARRAY(data) + offset, datatype, pageSize);
    }

    void
    getPageFromT(size_t offset, char * page, DataType datatype, size_t pageSize)
    {
        castPage2Datatype(MULTIDIM_ARRAY(data) + offset, page, datatype,
                          pageSize);
    }

    void
    getCastConvertPageFromT(size_t offset, char * page, DataType datatype,
                            size_t pageSize, double min0, double max0, CastWriteMode castMode =
    CW_CONVERT) const
    {
        castConvertPage2Datatype(MULTIDIM_ARRAY(data) + offset, page, datatype,
                                 pageSize, min0, max0, castMode);
    }

    /** Check if file Datatype is the same as the declared image object (T type) to use mmap.
     */
    bool
    checkMmapT(DataType datatype) override
    {

        switch (datatype)
        {
            case DT_Unknown:
                REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR: datatype is Unknown_Type");
            case DT_UHalfByte:
                return 0;
            case DT_UChar:
                return typeid(T) == typeid(unsigned char);
            case DT_SChar:
                return typeid(T) == typeid(char);
            case DT_UShort:
                return typeid(T) == typeid(unsigned short);
            case DT_Short:
                return typeid(T) == typeid(short);
            case DT_UInt:
                return typeid(T) == typeid(unsigned int);
            case DT_Int:
                return typeid(T) == typeid(int);
            case DT_Long:
                return typeid(T) == typeid(long);
            case DT_Float:
                return typeid(T) == typeid(float);
            case DT_Double:
                return typeid(T) == typeid(double);
            case DT_HalfFloat:
                return typeid(T) == typeid(half_float::half);
            default:
            {
                std::cerr << "Datatype= " << datatype << std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT, " ERROR: cannot cast datatype to T");
            }
        }
    }

    /** flip image around X axis
     *
     */
    void
    mirrorY(void)
    {
        T aux(0);
        size_t Z, Y, X, N, Y2;

        X = XSIZE(data);
        Y = YSIZE(data);
        Z = ZSIZE(data);
        N = NSIZE(data);
        Y2 = Y / 2;
        Y--;
        for (size_t l = 0; l < N; ++l)
            for (size_t k = 0; k < Z; ++k)
                for (size_t i = 0; i < Y2; ++i)
                    for (size_t j = 0; j < X; ++j)
                    {
                        aux = DIRECT_NZYX_ELEM(data, l, k, i, j);
                        DIRECT_NZYX_ELEM(data, l, k, i, j) = DIRECT_NZYX_ELEM(data, l, k,
                                                                              Y - i, j);
                        DIRECT_NZYX_ELEM(data, l, k, Y-i, j) = aux;
                    }
    }

    /** flip image around Y axis
     *
     */
    void
    mirrorX(void)
    {
        T aux(0);
        size_t Z, Y, X, N, X2;

        X = XSIZE(data);
        Y = YSIZE(data);
        Z = ZSIZE(data);
        N = NSIZE(data);
        X2 = X / 2;
        X--;
        for (size_t l = 0; l < N; ++l)
            for (size_t k = 0; k < Z; ++k)
                for (size_t i = 0; i < Y; ++i)
                    for (size_t j = 0; j < X2; ++j)
                    {
                        aux = DIRECT_NZYX_ELEM(data, l, k, i, j);
                        DIRECT_NZYX_ELEM(data, l, k, i, j) = DIRECT_NZYX_ELEM(data, l, k,
                                                                              i, X - j);
                        DIRECT_NZYX_ELEM(data, l, k, i, X - j) = aux;
                    }
    }

    void
    selfApplyGeometry(int SplineDegree, bool wrap = xmipp_transformation::WRAP,
                      bool only_apply_shifts = false);

    /* Read an image with a lower resolution as a preview image.
     * If Zdim parameter is not passed, then all slices are rescaled.
     * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
     */
    int
    readPreview(const FileName &name, size_t Xdim, size_t Ydim = 0,
                int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);

    /** Returns an image with a lower resolution as a preview image.
     * If Zdim parameter is not passed, then all slices are rescaled.
     * If Ydim is not passed, then Ydim is rescaled same factor as Xdim.
     */
    void
    getPreview(ImageBase *imgBOut, size_t Xdim, size_t Ydim = 0,
               int select_slice = CENTRAL_SLICE, size_t select_img = FIRST_IMAGE);

    /** It changes the behavior of the internal multidimarray so it points to a specific slice/image
     *  from a stack, volume or stack of volumes. No information is deallocated from memory, so it is
     *  also possible to repoint to the whole stack,volume... (passing select_slice = ALL_SLICES and
     *  selec_img = ALL_IMAGES).
     *
     *  The options for select_slice are:
     *
     *    - a slice number,
     *    - CENTRAL_SLICE, to automatically select the central slice of the volume,
     *    - ALL_SLICES, to recover the whole volume.
     *
     *  The options for selec_img are:
     *
     *    - a image number of the stack,
     *    - ALL_IMAGES, to recover the whole stack.
     *
     *  If a specific slice number is selected, then a specific image from the stack must be
     *  also selected. Otherwise, FIRST_IMAGE is proposed.
     *
     *  If Image Object is read using readPreview method, movePointerTo only works when rescaling
     *  the image in X-Y plane only, but all slices must be read.
     */
    void
    movePointerTo(int select_slice = ALL_SLICES, size_t select_img = ALL_IMAGES)
    {
        if (MULTIDIM_ARRAY(VOLMATRIX(*this)) == NULL)
            REPORT_ERROR(ERR_MULTIDIM_EMPTY,
                         "Image::movePointerTo: Image is empty");
        if (select_slice > (int) aDimFile.zdim)
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         formatString(
                                 "movePointerTo: Selected slice %4d cannot be higher than Z size %4d.",
                                 select_slice, aDimFile.zdim));
        else if (select_img > aDimFile.ndim)
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         formatString(
                                 "movePointerTo: Selected image %4d cannot be higher than N size %4d.",
                                 select_img, aDimFile.ndim));

        ArrayDim newDim;
        VOLMATRIX(*this).getDimensions(newDim);

        /* Restore the real dimensions from the MDA, as it may be
         * X-Y dimensioned different from file (readPreview). */
        newDim.zdim = aDimFile.zdim;
        newDim.ndim = aDimFile.ndim;

        int phys_slice;

        switch (select_slice)
        {
            case CENTRAL_SLICE:
                phys_slice = aDimFile.zdim / 2;
                newDim.zdim = 1;
                break;
            case ALL_SLICES:
                phys_slice = 0;
                break;
            default:
                phys_slice = select_slice - 1;
                newDim.zdim = 1;
                break;
        }

        /* If we select a single slice, we are forced to chose also an image.
         * as at this moment we cannot select the same slice from different images at a time.
         */
        if (select_slice > 0 && select_img == ALL_IMAGES)
            select_img = FIRST_IMAGE;
        if (select_img > ALL_IMAGES)
            newDim.ndim = 1;

        VOLMATRIX(*this).setDimensions(newDim);

        size_t newVirtualOffset = YXSIZE(VOLMATRIX(*this))
                                  * (aDimFile.zdim * IMG_INDEX(select_img) + phys_slice);
        MULTIDIM_ARRAY(VOLMATRIX(*this)) += (newVirtualOffset - virtualOffset);
        virtualOffset = newVirtualOffset;
    }

    /** Write an entire page as datatype
     *
     * A page of datasize_n elements T is cast to datatype and written to fimg
     * The memory for the casted page is allocated and freed internally.
     */
    void
    writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n)
    {
        size_t datasize = datasize_n * gettypesize(datatype);
        char * fdata = (char *) askMemory(datasize);
        castPage2Datatype(MULTIDIM_ARRAY(data), fdata, datatype, datasize_n);
        fwrite(fdata, datasize, 1, fimg);
        freeMemory(fdata, datasize);
    }

    /** Copy the MDA and the fields related to the (possible) original file
     */
    Image<T>&
    operator=(const Image<T> &op1)
    {
        this->copy(op1);
        data = op1.data;
        return *this;
    }

    /** Data access
     *
     * This operator can be used to access the data multidimarray.
     * In this way we could resize an image just by
     * resizing its associated matrix or we could add two images by adding their
     * matrices.
     * @code
     * I().resize(128, 128);
     * I2() = I1() + I2();
     * @endcode
     */
    MultidimArray<T>&
    operator()()
    {
        return data;
    }
    const MultidimArray<T>&
    operator()() const
    {
        return data;
    }

    /** Pixel access
     *
     * This operator is used to access a pixel within a 2D image. This is a
     * logical access, so you could access to negative positions if the image
     * has been defined so (see the general explanation for the class).
     *
     * @code
     * std::cout << "Grey level of pixel (-3,-3) of the image = " << I(-3, -3)
     * << std::endl;
     *
     * I(-3, -3) = I(-3, -2);
     * @endcode
     */
    T&
    operator()(int i, int j) const
    {
        return A2D_ELEM(data, i, j);
    }
    /**
     * equal operator
     */
    bool
    operator==(const Image<T> &i1) const
    {
        return (this->data == i1.data);
    }

    /** Voxel access
     *
     * This operator is used to access a voxel within a 3D image. This is a
     * logical access, so you could access to negative positions if the image
     * has been defined so (see the general explanation for the class).
     *
     * @code
     * std::cout << "Grey level of pixel (-3,-3, 1) of the volume = " << I(-3, -3, 1)
     * << std::endl;
     *
     * I(-3, -3, 1) = I(-3, -2, 0);
     * @endcode
     */
    T&
    operator()(int k, int i, int j) const
    {
        return A3D_ELEM(data, k, i, j);
    }

    /** Get Image dimensions
     */
    //    void getDimensions(int &Xdim, int &Ydim, int &Zdim, size_t &Ndim) const
    //    {
    //        Xdim = XSIZE(data);
    //        Ydim = YSIZE(data);
    //        Zdim = ZSIZE(data);
    //        Ndim = NSIZE(data);
    //    }
    size_t
    getSize() const
    {
        return NZYXSIZE(data);
    }

    /** Get geometric transformation matrix from 2D-image header
     */
    void
    getTransformationMatrix(Matrix2D<double> &A, bool only_apply_shifts = false,
                            const size_t n = 0);

    /** Sum this object with other file and keep in this object
     */
    void
    sumWithFile(const FileName &fn)
    {
        Image<T> aux;
        aux.read(fn, DATA, ALL_IMAGES, true);
        (*this)() += aux();
    }

    /**
     *  Specific read functions for different file formats
     */
    //#include "rwTIFF.h"
protected:

    /** Apply geometry in referring metadata to the image */
    void applyGeo(const MDRow &row, bool only_apply_shifts = false, bool wrap = xmipp_transformation::WRAP) override;

    //apply geo has not been defined for volumes
    //and only make sense when reading data

    /** Set the image dimensions
     */
    void
    setDimensions(int Xdim, int Ydim, int Zdim, size_t Ndim)
    {
        data.setDimensions(Xdim, Ydim, Zdim, Ndim);
        data.getDimensions(aDimFile);
    }

private:

    void setDimensions(ArrayDim &aDim) override
    {
        ImageBase::setDimensions(aDim);
    }

    static bool isValidAxisOrder(const std::array<int, 4>& order)
    {
        std::set<int> uniqueValues;
        
        for (int value : order) {
            // Check if the value is not in the range [0, 3] or is not unique.
            if (value < 0 || value > 3 || !uniqueValues.insert(value).second)
                return false;
        }

        return true;
    }

    /** Read the raw data
     */
    void
    readData(FILE* fimg, size_t select_img, DataType datatype, size_t pad) override
    {
        //#define DEBUG
#ifdef DEBUG
        std::cerr<<"entering readdata"<<std::endl;
      std::cerr<<" readData flag= "<<dataMode<<std::endl;
#endif
#undef DEBUG


        if(!isValidAxisOrder(axisOrder))
        {
            reportWarning("Image::readData: Invalid axis ordering. Defaulting to 0,1,2,3 ");
            axisOrder = defaultAxisOrder;
        }

        if (dataMode < DATA)
        {
            if (axisOrder != defaultAxisOrder) {
                std::array<size_t, 4> sizes;
                getDimensions(sizes[3], sizes[2], sizes[1], sizes[0]);
                transposeAxisSizes(sizes, axisOrder, sizes);
                setDimensions(sizes[3], sizes[2], sizes[1], sizes[0]);
            } 
            return;
        }

        if (datatype == DT_UHalfByte){
            //REPORT_ERROR(ERR_MMAP, "Image Class::readData not supported for  "
            //                       "data type " + datatype2Str(DT_UHalfByte));
            std::cout<<"redirecting from readData to readData4bits!"<<std::endl;
            readData4bit(fimg, select_img, datatype, pad);
        }
        // If only half of a transform is stored, it needs to be handled
        if (transform == Hermitian || transform == CentHerm)
            data.setXdim(XSIZE(data) / 2 + 1);

        size_t selectImgOffset, readsize, readsize_n, pagemax = 4194304; //4Mb
        size_t datatypesize = gettypesize(datatype);
        size_t pagesize = ZYXSIZE(data) * datatypesize;
        size_t haveread_n = 0;

        selectImgOffset = offset + IMG_INDEX(select_img) * (pagesize + pad);

        // Flag to know that data is not going to be mapped although mmapOn is true
        if (mmapOnRead && (!checkMmapT(datatype) || swap > 0 || axisOrder != defaultAxisOrder))
        {
            String warnMessage;
            if (swap > 0)
                reportWarning("Image::readData: File endianness is swapped and not "
                              "compatible with mmap. Loading into memory.");
            else if (axisOrder != defaultAxisOrder)
                reportWarning("Image::readData: Axis order is not standard 0,1,2,3, which makes it "
                                "incompatible with memory mapping. Loading into memory.");
            else
                reportWarning(
                        "Image::readData: File datatype and image declaration not "
                        "compatible with mmap. Loading into memory.");

            mmapOnRead = false;
            mFd = -1;
        }

        if (mmapOnRead)
        {
            // Image mmapOn is not compatible with Multidimarray mmapOn
            if (data.mmapOn)
                REPORT_ERROR(ERR_MULTIDIM_DIM,
                             "Image Class::ReadData: mmap option can not be selected simultaneously\
                             for both Image class and its Multidimarray.");
            if ( NSIZE(data) > 1)
            {
                REPORT_ERROR(ERR_MMAP, "Image Class::ReadData: mmap with multiple "
                                       "images file not compatible. Try selecting a unique image.");
            }
            mappedOffset = selectImgOffset;
            mappedSize = mappedOffset + pagesize;
            mmapFile();
        }
        else
        {
            // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
            //if memory already allocated use it (no resize allowed)
            data.coreAllocateReuse();
            //ROB
// #define DEBUG
#ifdef DEBUG

            data.printShape();
        printf("DEBUG: Page size: %ld offset= %ld \n", pagesize, offset);
        printf("DEBUG: Swap = %d  Pad = %ld  Offset = %ld\n", swap, pad, offset);
        printf("DEBUG: myoffset = %ld select_img= %ld \n", selectImgOffset, select_img);
        printf("DEBUG: NSIZE = %lu \n", NSIZE(data));
#endif
#undef DEBUG

        if (checkMmapT(datatype) && !swap && axisOrder == defaultAxisOrder) {
            // printf( "type is same, reading without cast\n" );

            size_t slice_elements = ZYXSIZE(data);

            if (fseek(fimg, selectImgOffset, SEEK_SET) == -1)
                REPORT_ERROR(ERR_IO_SIZE, "readData: can not seek the file pointer");

            if (pad == 0) {
                // printf( "pad == 0, reading with one fread \n" );
                if (fread(MULTIDIM_ARRAY(data), pagesize * NSIZE(data), 1, fimg) != 1) {
                        REPORT_ERROR(ERR_IO_NOREAD, "readData: cannot read the whole image slice");
                    }
            } else {
                for (size_t n = 0; n < NSIZE(data); ++n) {
                    if (fread(MULTIDIM_ARRAY(data) + slice_elements * n, pagesize, 1, fimg) != 1) {
                        REPORT_ERROR(ERR_IO_NOREAD, "readData: cannot read the whole image slice");
                    }

                    if (fseek(fimg, pad, SEEK_CUR) == -1) {
                            REPORT_ERROR(ERR_IO_SIZE, "readData: cannot seek the file pointer");
                    }
                }
            }

            
        } else {
            // std::cout << "original read" << std::endl;
            char* page = NULL;

            if (pagesize > pagemax)
                page = (char *) askMemory(pagemax * sizeof(char));
            else
                page = (char *) askMemory(pagesize * sizeof(char));

            if (fseek(fimg, selectImgOffset, SEEK_SET) == -1)
                REPORT_ERROR(ERR_IO_SIZE, "readData: can not seek the file pointer");
            for (size_t myn = 0; myn < NSIZE(data); myn++)
            {
                for (size_t myj = 0; myj < pagesize; myj += pagemax) //pagesize size of object
                {
                    // Read next page. Divide pages larger than pagemax
                    readsize = pagesize - myj;
                    if (readsize > pagemax)
                        readsize = pagemax;
                    readsize_n = readsize / datatypesize;

                    //Read page from disc
                    if (fread(page, readsize, 1, fimg) != 1)
                        REPORT_ERROR(ERR_IO_NOREAD, "Cannot read the whole page");
                    //swap per page
                    if (swap)
                        swapPage(page, readsize, datatype, swap);
                    // cast to T per page
                    castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, datatype,
                               readsize_n);
                    haveread_n += readsize_n;
                }
                if (pad > 0)
                    //fread( padpage, pad, 1, fimg);
                    if (fseek(fimg, pad, SEEK_CUR) == -1)
                        REPORT_ERROR(ERR_IO_SIZE,
                                     "readData: can not seek the file pointer");
            }
            // Transposing multidim array
            if (axisOrder != defaultAxisOrder) {
                transposeInPlace(data, axisOrder);
                data.getDimensions(aDimFile);
            }
            //if ( pad > 0 )
            //    freeMemory(padpage, pad*sizeof(char));
            if (page)
                freeMemory(page, pagesize * sizeof(char));

        }

#ifdef DEBUG

            printf("DEBUG img_read_data: Finished reading and converting data\n");
#endif

        }
        return;
    }



    /** Read the raw data from compressed 4bit images
     * We are assuming the values are stored in 4bits (2 values in 1 byte)
    */
    void
    readData4bit(FILE* fimg, size_t select_img, DataType datatype, size_t pad)
    {

        if (dataMode < DATA)
            return;

        if (datatype != DT_UHalfByte){
            REPORT_ERROR(ERR_MMAP, "Image Class::readData4bit  not supported for  "
                                   "data type different than " + datatype2Str(DT_UHalfByte));
        }

        size_t selectImgOffset; //4Mb
        size_t itemSize = ZYXSIZE(data);
        size_t pagesizeF = itemSize /2;
        size_t pagesizeM = itemSize;
        //size_t pagesizeHalf = pagesize/2;
        size_t haveread_n = 0;

        char* page = NULL;    // Compressed

        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)
        data.coreAllocateReuse();

        page = (char *) askMemory(pagesizeM);

        // Calculate the staring reading pointer.
        selectImgOffset = offset + IMG_INDEX(select_img) * (pagesizeF + pad);

        if (fseek(fimg, selectImgOffset, SEEK_SET) == -1)
            REPORT_ERROR(ERR_IO_SIZE, "readData4bit: can not seek the file pointer");
        for (size_t myn = 0; myn < NSIZE(data); myn++)
        {

            //Read page from disc
            if (fread(page+pagesizeF, pagesizeF, 1, fimg) != 1)
                REPORT_ERROR(ERR_IO_NOREAD, "Cannot read the whole page");

            // cast to T per page

            size_t start = pagesizeF;
            uint8_t mask = 15; // 00001111

            for (size_t i = 0, j = start; i < pagesizeM - 1; i += 2, ++j)
            {
                char& value = *(page+j);
                page[i] = value & mask; // take the lower 4 bits
                page[i+1] = (value >> 4) & mask; // take the upper 4 bits
            }

            castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, datatype,
                       pagesizeM);
            haveread_n += pagesizeM;

            if (pad > 0)
                //fread( padpage, pad, 1, fimg);
                if (fseek(fimg, pad, SEEK_CUR) == -1)
                    REPORT_ERROR(ERR_IO_SIZE,
                                 "readData4bit: can not seek the file pointer");

            //if ( pad > 0 )
            //    freeMemory(padpage, pad*sizeof(char));

#ifdef DEBUG

            printf("DEBUG readData4bit: Finished reading and converting data\n");
#endif
        }
        if (page)
            freeMemory(page, pagesizeM * sizeof(char));

        return;
    }








    /* Write the raw date after a data type casting.
     */
    void
    writeData(FILE* fimg, size_t offset, DataType wDType, size_t datasize_n,
              CastWriteMode castMode = CW_CAST)
    {
        size_t dTypeSize = gettypesize(wDType);
        size_t datasize = datasize_n * dTypeSize;
        size_t ds2Write = rw_max_page_size;
        size_t dsN2Write = rw_max_page_size / dTypeSize;
        size_t rw_max_n = dsN2Write;

        char* fdata;
        double min0 = 0, max0 = 0;

        if (wDType == myT() && castMode == CW_CONVERT)
            castMode = CW_CAST;

        if (castMode != CW_CAST)
            data.computeDoubleMinMaxRange(min0, max0, offset, datasize_n);

        if ( checkMmapT(wDType) ) {
            fwrite( MULTIDIM_ARRAY(data) + offset, datasize, 1, fimg );
            return;
        }

        if (datasize > rw_max_page_size)
            fdata = (char *) askMemory(rw_max_page_size * sizeof(char));
        else
            fdata = (char *) askMemory(datasize * sizeof(char));

        for (size_t writtenDataN = 0; writtenDataN < datasize_n; writtenDataN +=
                                                                         rw_max_n)
        {

            if (writtenDataN + rw_max_n > datasize_n)
            {
                dsN2Write = datasize_n - writtenDataN;
                ds2Write = dsN2Write * dTypeSize;
            }

            if (castMode == CW_CAST)
                castPage2Datatype(MULTIDIM_ARRAY(data) + offset + writtenDataN, fdata,
                                  wDType, dsN2Write);
            else
                castConvertPage2Datatype(MULTIDIM_ARRAY(data) + offset + writtenDataN,
                                         fdata, wDType, dsN2Write, min0, max0, castMode);

            //swap per page
            if (swapWrite)
                swapPage(fdata, ds2Write, wDType);

            fwrite(fdata, ds2Write, 1, fimg);
        }
        freeMemory(fdata, rw_max_page_size);
    }

    /* Mmap the Image class to an image file.
     */
    void
    mmapFile();

    /* Munmap the image file.
     */
    void
    munmapFile()
    {
#ifdef XMIPP_MMAP
        munmap((char*) (data.data) - mappedOffset, mappedSize);
        close(mFd);
        data.data = NULL;
        mappedSize = mappedOffset = 0;
#else

        REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

    }

    /* Return the datatype of the current image object
     */
    DataType
    myT() const
    {
        if (typeid(T) == typeid(unsigned char))
            return DT_UChar;
        else if (typeid(T) == typeid(char))
            return DT_SChar;
        else if (typeid(T) == typeid(unsigned short))
            return DT_UShort;
        else if (typeid(T) == typeid(short))
            return DT_Short;
        else if (typeid(T) == typeid(unsigned int))
            return DT_UInt;
        else if (typeid(T) == typeid(int))
            return DT_Int;
        else if (typeid(T) == typeid(unsigned int))
            return DT_UInt;
        else if (typeid(T) == typeid(int))
            return DT_Int;
        else if (typeid(T) == typeid(long))
            return DT_Long;
        else if (typeid(T) == typeid(float))
            return DT_Float;
        else if (typeid(T) == typeid(double))
            return DT_Double;
        else if (typeid(T) == typeid(std::complex<short>))
            return DT_CShort;
        else if (typeid(T) == typeid(std::complex<int>))
            return DT_CInt;
        else if (typeid(T) == typeid(std::complex<float>))
            return DT_CFloat;
        else if (typeid(T) == typeid(std::complex<double>))
            return DT_CDouble;
        else if (typeid(T) == typeid(bool))
            return DT_HalfFloat;
        else if (typeid(T) == typeid(half_float::half))
            return DT_Bool;
        else
            return DT_Unknown;
    }

    /* friend declaration for stacks handling purposes
     */
    friend class ImageCollection;
    template<typename TT>
    friend class Image;

}
;

// Special cases for complex numbers
template<>
void
Image<std::complex<double> >::castPage2T(char * page,
                                         std::complex<double> * ptrDest, DataType datatype, size_t pageSize);
template<>
void
Image<std::complex<double> >::castPage2Datatype(std::complex<double> * srcPtr,
                                                char * page, DataType datatype, size_t pageSize) const;
template<>
void
Image<std::complex<double> >::castConvertPage2Datatype(
        std::complex<double> * srcPtr, char * page, DataType datatype,
        size_t pageSize, double min0, double max0, CastWriteMode castMode) const;

//@}
#endif
