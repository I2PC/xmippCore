/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include "xmipp_image.h"
#include "xmipp_image_generic.h"
#include "matrix2d.h"
#include "transformations.h"
#include <sys/stat.h>
#include "metadata_static.h"

template<typename T>
int Image<T>::readPreview(const FileName &name, size_t Xdim, size_t Ydim,
                int select_slice, size_t select_img) // FIXME this should be moved to image_generic.*
{
    // Zdim is used to choose the slices: -1 = CENTRAL_SLICE, 0 = ALL_SLICES, else This Slice

    ImageGeneric im;
    size_t imXdim, imYdim, imZdim, Zdim;
    int err;
    err = im.readMapped(name, select_img);
    im.getDimensions(imXdim, imYdim, imZdim);
    ImageInfo imgInfo;
    im.getInfo(imgInfo);

    //Set information from image file
    setName(name);
    setDatatype(imgInfo.datatype);
    aDimFile = imgInfo.adim;

    im().setXmippOrigin();

    double scale;

    // If only Xdim is passed, it is the higher allowable size, for any dimension
    if (Ydim == 0 && imXdim < imYdim)
    {
        Ydim = Xdim;
        scale = ((double) Ydim) / ((double) imYdim);
        Xdim = (int) (scale * imXdim);
    }
    else
    {
        scale = ((double) Xdim) / ((double) imXdim);
        if (Ydim == 0)
            Ydim = (int) (scale * imYdim);
    }

    int mode = (scale <= 1) ? xmipp_transformation::NEAREST : xmipp_transformation::LINEAR; // If scale factor is higher than 1, LINEAR mode is used to avoid artifacts

    if (select_slice > ALL_SLICES) // In this case a specific slice number has been chosen (Not central slice)
    {
        MultidimArrayGeneric array(im(), select_slice - 1);
        array.setXmippOrigin();

        scaleToSize(mode, IMGMATRIX(*this), array, Xdim, Ydim);
    }
    else // Otherwise, All slices or Central slice is selected
    {
        Zdim = (select_slice == ALL_SLICES) ? imZdim : 1;
        scaleToSize(mode, IMGMATRIX(*this), im(), Xdim, Ydim, Zdim);
    }

    IMGMATRIX(*this).resetOrigin();
    return err;
}

template<typename T>
void Image<T>::mmapFile()
    {
#ifdef XMIPP_MMAP
        if (this->hFile->mode == WRITE_READONLY)
            mFd = open(dataFName.c_str(), O_RDONLY, S_IREAD);
        else
            mFd = open(dataFName.c_str(), O_RDWR, S_IREAD | S_IWRITE);

        if (mFd == -1)
        {
            if (errno == EACCES)
                REPORT_ERROR(ERR_IO_NOPERM,
                             formatString(
                                     "Image Class::mmapFile: permission denied when opening %s",
                                     dataFName.c_str()));
            else
                REPORT_ERROR(ERR_IO_NOTOPEN,
                             "Image Class::mmapFile: Error opening the image file to be mapped.");
        }
        char * map;
        const size_t pagesize = sysconf(_SC_PAGESIZE);
        size_t offsetPages = (mappedOffset / pagesize) * pagesize;
        mappedOffset -= offsetPages;
        mappedSize -= offsetPages;

        if (this->hFile->mode == WRITE_READONLY)
            map = (char*) mmap(0, mappedSize, PROT_READ, MAP_SHARED, mFd,
                               offsetPages);
        else
            map = (char*) mmap(0, mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED,
                               mFd, offsetPages);

        if (map == MAP_FAILED)
            REPORT_ERROR(ERR_MMAP_NOTADDR,
                         formatString("Image Class::mmapFile: mmap of image file failed. Error: %s", strerror(errno)));
        data.data = reinterpret_cast<T*>(map + mappedOffset);
        data.nzyxdimAlloc = XSIZE(data) * YSIZE(data) * ZSIZE(data) * NSIZE(data);
#else

        REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

    }


// Special cases for complex numbers
template<>
void Image< std::complex< double > >::castPage2T(char * page,
        std::complex<double> * ptrDest,
        DataType datatype,
        size_t pageSize)
{

    switch (datatype)
    {
    case DT_CShort:
        {
            std::complex<short> * ptr = (std::complex<short> *) page;
            for(size_t i=0; i<pageSize;i++)
                ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
        }
        break;
    case DT_CInt:
		{
			std::complex<int> * ptr = (std::complex<int> *) page;
			for(size_t i=0; i<pageSize;i++)
			ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
		}
		break;
    case DT_CFloat:
        {
			std::complex<float> * ptr = (std::complex<float> *) page;
			for(size_t i=0; i<pageSize;i++)
				ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
		}
		break;
    case DT_CDouble:
            memcpy(ptrDest, page, pageSize*sizeof(std::complex<double>));
        break;
    default:
		std::cerr<<"Datatype= "<<datatype<<std::endl;
		REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to std::complex<double>");
        break;
    }
}

template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex<double> * srcPtr,
        char * page,
        DataType datatype,
        size_t pageSize) const
{
    switch (datatype)
    {
    case DT_CShort:
        {
            short  * ptr = (short *) page;
            double * srcPtrd = (double *)srcPtr;
            for(size_t i=0; i<pageSize;i++)
            {
                *ptr=(short) *srcPtrd; ptr++; srcPtrd++;
                *ptr=(short) *srcPtrd; ptr++; srcPtrd++;
            }
        }
        break;
    case DT_CInt:
		{
                        int  * ptr = (int *) page;
                        double * srcPtrd = (double *)srcPtr;
			for(size_t i=0; i<pageSize;i++)
                        {
                            *ptr=(int) *srcPtrd; ptr++; srcPtrd++;
                            *ptr=(int) *srcPtrd; ptr++; srcPtrd++;
                        }
		}
		break;
	case DT_CFloat:
		{
			std::complex<float> * ptr = (std::complex<float> *) page;
			for(size_t i=0; i<pageSize;i++)
				ptr[i] = (std::complex<float>)srcPtr[i];
		}
		break;
	case DT_CDouble:
		memcpy(page, srcPtr, pageSize*sizeof(std::complex<double>));
		break;
	default:
		REPORT_ERROR(ERR_TYPE_INCORRECT,formatString("ERROR: cannot cast type number %d to complex<double>",datatype));
		break;
	}
}

template<>
void Image< std::complex< double > >::castConvertPage2Datatype(std::complex< double > * srcPtr,
        char * page, DataType datatype, size_t pageSize,double min0,double max0,CastWriteMode castMode) const
{

    switch (datatype)
    {
    case DT_CFloat:
        {
            std::complex<float> * ptr = (std::complex<float> *) page;
            for(size_t i=0; i<pageSize;i++)
                ptr[i] = (std::complex<float>)srcPtr[i];
        }
        break;
    default:
		REPORT_ERROR(ERR_TYPE_INCORRECT,formatString("ERROR: cannot cast&convert type number %d to complex<double>",datatype));
		break;
    }
}

template<typename T>
void Image<T>::selfApplyGeometry(int SplineDegree, bool wrap,
                  bool only_apply_shifts)
{
    //apply geo has not been defined for volumes
    //and only make sense when reading data
    if (data.getDim() < 3 && dataMode >= DATA)
    {
        Matrix2D<double> A;
        getTransformationMatrix(A, only_apply_shifts);
        if (!A.isIdentity())
        {
            MultidimArray<T> tmp = MULTIDIM_ARRAY(*this);
            applyGeometry(SplineDegree, MULTIDIM_ARRAY(*this), tmp, A, xmipp_transformation::IS_NOT_INV,
                          wrap);
        }
    }
}

template<typename T>
void Image<T>::getTransformationMatrix(Matrix2D<double> &A, bool only_apply_shifts,
                        const size_t n)
{
    // This has only been implemented for 2D images...
    MULTIDIM_ARRAY(*this).checkDimension(2);
    A.resizeNoCopy(3, 3);
    geo2TransformationMatrix(*MD[n], A, only_apply_shifts);
}

template<typename T>
void Image<T>::getPreview(ImageBase *imgBOut, size_t Xdim, size_t Ydim,
           int select_slice, size_t select_img)
{
    // Zdim is used to choose the slices: -1 = CENTRAL_SLICE, 0 = ALL_SLICES, else This Slice

    size_t Zdim;
    ArrayDim imAdim;
    MULTIDIM_ARRAY(*this).getDimensions(imAdim);
    MULTIDIM_ARRAY(*this).setXmippOrigin();

    double scale;

    // If only Xdim is passed, it is the higher allowable size, for any dimension
    if (Ydim == 0 && imAdim.xdim < imAdim.ydim)
    {
        Ydim = Xdim;
        scale = ((double) Ydim) / ((double) imAdim.ydim);
        Xdim = (int) (scale * imAdim.xdim);
    }
    else
    {
        scale = ((double) Xdim) / ((double) imAdim.xdim);
        if (Ydim == 0)
            Ydim = (int) (scale * imAdim.ydim);
    }

    Image<T> &imgOut = *((Image<T>*) imgBOut);

    int mode = (scale <= 1) ? xmipp_transformation::NEAREST : xmipp_transformation::LINEAR; // If scale factor is higher than 1, LINEAR mode is used to avoid artifacts

    if (select_slice > ALL_SLICES) // In this case a specific slice number has been chosen (Not central slice)
    {
        movePointerTo(select_slice, select_img);
        scaleToSize(mode, IMGMATRIX(imgOut), IMGMATRIX(*this), Xdim, Ydim);
    }
    else // Otherwise, All slices or Central slice is selected
    {
        movePointerTo(ALL_SLICES, select_img);
        Zdim = (select_slice == ALL_SLICES) ? imAdim.zdim : 1;
        scaleToSize(mode, IMGMATRIX(imgOut), IMGMATRIX(*this), Xdim, Ydim,
                    Zdim);
    }

    movePointerTo();
    IMGMATRIX(*this).resetOrigin();

    // We set the actual dimesions of th MDA to the imageOut as if it were read from file.
    imgOut.setADimFile(IMGMATRIX(imgOut).getDimensions());
}

template<typename T>
void Image<T>::applyGeo(const MDRow &row, bool only_apply_shifts, bool wrap) {
    //This implementation does not handle stacks,
    //read in a block
    if (data.ndim != 1)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     "Geometric transformation cannot be applied to stacks!!!");

    if (MD.size() == 0)
        MD.push_back(std::unique_ptr<MDRowVec>(new MDRowVec(MDL::emptyHeaderVec())));
    MDRow &rowAux = *MD[0];

    if (!row.containsLabel(MDL_TRANSFORM_MATRIX))
    {
        double aux;
        //origins
        if (row.getValue(MDL_ORIGIN_X, aux))
            rowAux.setValue(MDL_ORIGIN_X, aux);
        if (row.getValue(MDL_ORIGIN_Y, aux))
            rowAux.setValue(MDL_ORIGIN_Y, aux);
        if (row.getValue(MDL_ORIGIN_Z, aux))
            rowAux.setValue(MDL_ORIGIN_Z, aux);
        //shifts
        if (row.getValue(MDL_SHIFT_X, aux))
            rowAux.setValue(MDL_SHIFT_X, aux);
        if (row.getValue(MDL_SHIFT_Y, aux))
            rowAux.setValue(MDL_SHIFT_Y, aux);
        if (row.getValue(MDL_SHIFT_Z, aux))
            rowAux.setValue(MDL_SHIFT_Z, aux);
        //rotations
        if (row.getValue(MDL_ANGLE_ROT, aux))
            rowAux.setValue(MDL_ANGLE_ROT, aux);
        if (row.getValue(MDL_ANGLE_TILT, aux))
            rowAux.setValue(MDL_ANGLE_TILT, aux);
        if (row.getValue(MDL_ANGLE_PSI, aux))
            rowAux.setValue(MDL_ANGLE_PSI, aux);
        //scale
        if (row.getValue(MDL_SCALE, aux))
            rowAux.setValue(MDL_SCALE, aux);
        //weight
        if (row.getValue(MDL_WEIGHT, aux))
            rowAux.setValue(MDL_WEIGHT, aux);
        bool auxBool;
        if (row.getValue(MDL_FLIP, auxBool))
            rowAux.setValue(MDL_FLIP, auxBool);
    }

    //apply geo has not been defined for volumes
    //and only make sense when reading data
    if (data.getDim() < 3 && dataMode >= DATA)
    {
        Matrix2D<double> A;
        if (!row.containsLabel(MDL_TRANSFORM_MATRIX))
            getTransformationMatrix(A, only_apply_shifts);
        else
        {
            String matrixStr;
            row.getValue(MDL_TRANSFORM_MATRIX, matrixStr);
            string2TransformationMatrix(matrixStr, A, 3);
        }

        if (!A.isIdentity())
        {
            MultidimArray<T> tmp = MULTIDIM_ARRAY(*this);
            applyGeometry(xmipp_transformation::BSPLINE3, MULTIDIM_ARRAY(*this), tmp, A, xmipp_transformation::IS_NOT_INV,
                          wrap);
        }
    }
}

//template int Image<std::complex<double> >::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<bool>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<int>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<short>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<float>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<unsigned int>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<double>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<char>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<unsigned short>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<unsigned char>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<unsigned long>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
//template int Image<long>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template class Image<std::complex<double> >;
template class Image<bool>;
template class Image<int>;
template class Image<short>;
template class Image<float>;
template class Image<unsigned int>;
template class Image<double>;
template class Image<char>;
template class Image<unsigned short>;
template class Image<unsigned char>;
template class Image<unsigned long>;
template class Image<long>;
template class Image<half_float::half>;
