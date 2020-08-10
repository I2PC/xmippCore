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

    int mode = (scale <= 1) ? NEAREST : LINEAR; // If scale factor is higher than 1, LINEAR mode is used to avoid artifacts

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

template int Image<std::complex<double> >::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<bool>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<int>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<short>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<float>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<unsigned int>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<double>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<char>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<unsigned short>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<unsigned char>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<unsigned long>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
template int Image<long>::readPreview(FileName const&, unsigned long, unsigned long, int, unsigned long);
