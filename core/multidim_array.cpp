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

#include "multidim_array.h"
#include "matrix2d.h"
#include "multidim_array_base.h"

template<typename T>
void MultidimArray<T>::getSliceAsMatrix(size_t k, Matrix2D<T> &m) const
{
    m.resizeNoCopy(YSIZE(*this),XSIZE(*this));
    memcpy(&MAT_ELEM(m,0,0),&A3D_ELEM(*this,k,0,0),YSIZE(*this),XSIZE(*this)*sizeof(double));
}

template<typename T>
MultidimArray<T>& MultidimArray<T>::operator=(const Matrix2D<T>& op1)
{
    resizeNoCopy(MAT_YSIZE(op1), MAT_XSIZE(op1));
    memcpy(data,MATRIX2D_ARRAY(op1), MAT_SIZE(op1)*sizeof(T));

    return *this;
}

template<typename T>
void MultidimArray<T>::copy(Matrix2D<T>& op1) const
{
    op1.resizeNoCopy(YSIZE(*this), XSIZE(*this));
    memcpy(MATRIX2D_ARRAY(op1), MULTIDIM_ARRAY(*this), MULTIDIM_SIZE(*this)*sizeof(T));
}

template<typename T>
void MultidimArray<T>::rangeAdjustLS(const MultidimArray<T> &example,
                 const MultidimArray<int> *mask)
{
    if (NZYXSIZE(*this) <= 0)
        return;

    Matrix2D<double> A;
    A.initZeros(2,2);
    Matrix1D<double> b;
    b.initZeros(2);

    size_t n;
    T *ptr=NULL;
    if (mask==NULL)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
			MAT_ELEM(A,0,0)+=1;
			MAT_ELEM(A,0,1)+=*ptr;
			MAT_ELEM(A,1,1)+=(*ptr)*(*ptr);

			VEC_ELEM(b,0)+=DIRECT_MULTIDIM_ELEM(example,n);
			VEC_ELEM(b,1)+=DIRECT_MULTIDIM_ELEM(example,n)*(*ptr);
        }
    }
    else
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (DIRECT_MULTIDIM_ELEM(*mask,n))
        {
			MAT_ELEM(A,0,0)+=1;
			MAT_ELEM(A,0,1)+=*ptr;
			MAT_ELEM(A,1,1)+=(*ptr)*(*ptr);

			VEC_ELEM(b,0)+=DIRECT_MULTIDIM_ELEM(example,n);
			VEC_ELEM(b,1)+=DIRECT_MULTIDIM_ELEM(example,n)*(*ptr);
        }
    }

    b=A.inv()*b;
	T a0=(T)VEC_ELEM(b,0);
	T a1=(T)VEC_ELEM(b,1);
	if (mask==NULL)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
		*ptr = a0+a1*(*ptr);
	else
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (DIRECT_MULTIDIM_ELEM(*mask,n))
			*ptr = a0+a1*(*ptr);
}


// Window in 2D -----------------------------------------------------------
template<typename T>
void window2D(const MultidimArray<T> &Ibig, MultidimArray<T> &Ismall,
        size_t y0, size_t x0, size_t yF, size_t xF)
{
    Ismall.resizeNoCopy(yF - y0 + 1, xF - x0 + 1);
    STARTINGY(Ismall) = y0;
    STARTINGX(Ismall) = x0;
    /*
	FOR_ALL_ELEMENTS_IN_ARRAY2D(Ismall)
		A2D_ELEM(Ismall, i, j) = A2D_ELEM(Ibig, i, j);
    */
    size_t sizeToCopy=XSIZE(Ismall)*sizeof(T);
    for (int y=y0; y<=yF; y++)
    	memcpy( &A2D_ELEM(Ismall,y,STARTINGX(Ismall)), &A2D_ELEM(Ibig,y,STARTINGX(Ismall)), sizeToCopy);
}

template void window2D(const MultidimArray<double> &Ibig, MultidimArray<double> &Ismall, size_t y0, size_t x0, size_t yF, size_t xF);
template void window2D(const MultidimArray<float> &Ibig, MultidimArray<float> &Ismall, size_t y0, size_t x0, size_t yF, size_t xF);

// Show a complex array ---------------------------------------------------
template<>
std::ostream& operator<<(std::ostream& ostrm,
                         const MultidimArray< std::complex<double> >& v)
{
    if (v.xdim == 0)
        ostrm << "NULL MultidimArray\n";
    else
        ostrm << std::endl;

    for (size_t l = 0; l < NSIZE(v); l++)
    {
        if (NSIZE(v)>1)
            ostrm << "Image No. " << l << std::endl;
        for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
        {
            if (ZSIZE(v)>1)
                ostrm << "Slice No. " << k << std::endl;
            for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
            {
                for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
                    ostrm << A3D_ELEM(v, k, i, j) << ' ';
                ostrm << std::endl;
            }
        }
    }

    return ostrm;
}

template<>
void MultidimArray< std::complex< double > >::computeDoubleMinMax(double& minval, double& maxval) const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::computeDoubleMinMax not implemented for complex.");
}
template<>
void MultidimArray< std::complex< double > >::computeDoubleMinMaxRange(double& minval, double& maxval, size_t pos, size_t size) const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::computeDoubleMinMax not implemented for complex.");
}
template<>
void MultidimArray< std::complex< double > >::rangeAdjust(std::complex< double > minF, std::complex< double > maxF)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::rangeAdjust not implemented for complex.");
}

template<>
double MultidimArray< std::complex< double > >::computeAvg() const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::computeAvg not implemented for complex.");
}

template<>
void MultidimArray< std::complex< double > >::maxIndex(size_t &lmax, int& kmax, int& imax, int& jmax) const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"MultidimArray::maxIndex not implemented for complex.");
}

// void MultidimArray<double>::selfNormalizeInterval(double minPerc, double maxPerc, int Npix)
// {
//     std::vector<double> randValues; // Vector with random chosen values

//     for(int i=0, i<Npix, i++)
//     {
//         size_t indx = (size_t)rnd_unif(0, MULTIDIM_SIZE(*this));
//         randValues.push_back(DIRECT_MULTIDIM_ELEM(*this,indx))
//     }
//     std::sort(randValues.begin(),randValues.end());

//     double m = randValues[(size_t)(minPerc*Npix)];
//     double M = randValues[(size_t)(minPerc*Npix)];

//     FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
//     *ptr = 2/(M-m)*(*ptr-m)-1;

// }

template<>
bool operator==(const MultidimArray< std::complex< double > >& op1, const MultidimArray< std::complex< double > >& op2)
{
    double accuracy = XMIPP_EQUAL_ACCURACY;
    if (! op1.sameShape(op2) || op1.data==NULL || op2.data == NULL)
        return false;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1)
    if (   fabs(DIRECT_MULTIDIM_ELEM(op1,n).real() -
                DIRECT_MULTIDIM_ELEM(op2,n).real() > accuracy)
           ||
           fabs(DIRECT_MULTIDIM_ELEM(op1,n).imag() -
                DIRECT_MULTIDIM_ELEM(op2,n).imag() > accuracy)
       )
        return false;
    return true;
}

template<>
void MultidimArray< std::complex< double > >::getReal(MultidimArray<double> & realImg) const
{
    if (NZYXSIZE(*this) == 0)
    {
        realImg.clear();
        return;
    }

    realImg.resizeNoCopy(*this);
    double * ptr1 = (double*) MULTIDIM_ARRAY(*this);

    // Unroll the loop
    const size_t unroll=4;
    size_t nmax=(NZYXSIZE(*this)/unroll)*unroll;
    for (size_t n=0; n<nmax; n+=unroll)
    {
        DIRECT_MULTIDIM_ELEM(realImg, n)   = static_cast<double>(*ptr1++);
        ptr1++;
        DIRECT_MULTIDIM_ELEM(realImg, n+1) = static_cast<double>(*(ptr1++));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(realImg, n+2) = static_cast<double>(*(ptr1++));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(realImg, n+3) = static_cast<double>(*(ptr1++));
        ptr1++;
    }
    // Do the remaining elements
    for (size_t n=nmax; n<NZYXSIZE(*this); ++n)
    {
        DIRECT_MULTIDIM_ELEM(realImg, n) = static_cast<double>(*ptr1++);
        ptr1++;
    }

}

template<>
void MultidimArray< std::complex< double > >::getImag(MultidimArray<double> & imagImg) const
{
    if (NZYXSIZE(*this) == 0)
    {
        imagImg.clear();
        return;
    }

    imagImg.resizeNoCopy(*this);
    double * ptr1 = (double*) MULTIDIM_ARRAY(*this);

    // Unroll the loop
    const size_t unroll=4;
    size_t nmax=(NZYXSIZE(*this)/unroll)*unroll;
    for (size_t n=0; n<nmax; n+=unroll)
    {
        DIRECT_MULTIDIM_ELEM(imagImg, n)   = static_cast<double>(*(++ptr1));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(imagImg, n+1) = static_cast<double>(*(++ptr1));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(imagImg, n+2) = static_cast<double>(*(++ptr1));
        ptr1++;
        DIRECT_MULTIDIM_ELEM(imagImg, n+3) = static_cast<double>(*(++ptr1));
        ptr1++;
    }
    // Do the remaining elements
    for (size_t n=nmax; n<NZYXSIZE(*this); ++n)
    {
        DIRECT_MULTIDIM_ELEM(imagImg, n) = static_cast<double>(*(++ptr1));
        ptr1++;
    }

}


template<>
double MultidimArray<double>::interpolatedElement2D(double x, double y, double outside_value) const
{
    double dx0 = floor(x);
    int x0=(int)dx0;
    double fx = x - dx0;
    int x1 = x0 + 1;
    double dy0 = floor(y);
    int y0=(int)dy0;
    double fy = y - dy0;
    int y1 = y0 + 1;

    int i0=STARTINGY(*this);
    int j0=STARTINGX(*this);
    int iF=FINISHINGY(*this);
    int jF=FINISHINGX(*this);

    /* Next code avoids checking two times some variables that are doubles.
     * The code before was:
        ASSIGNVAL2D(d00,y0,x0);
        ASSIGNVAL2D(d01,y0,x1);
        ASSIGNVAL2D(d10,y1,x0);
        ASSIGNVAL2D(d11,y1,x1);
    */
    double d00, d10, d11, d01, *ref;
    if ((x0 >= j0) && (x1 <= jF) && (y0 >= i0) && (y1 <= iF))
    {
    	ref = &A2D_ELEM(*this, y0, x0);
		d00 = (*ref++);
    	d01 = (*ref);
    	ref = &A2D_ELEM(*this, y1, x0);
		d10 = (*ref++);
		d11 = (*ref);
    }
    else
    {
		bool outX0,outX1,outY0,outY1;
		outX0 = (x0 < j0) || (x0 > jF);
		outX1 = (x1 < j0) || (x1 > jF);
		outY0 = (y0 < i0) || (y0 > iF);
		outY1 = (y1 < i0) || (y1 > iF);

		if (outY0)
		{
			d00 = outside_value;
			d01 = outside_value;
			d10 = outside_value;
			d11 = outside_value;
		}
		else if(outY1)
		{
			d10 = outside_value;
			d11 = outside_value;

			if (outX0)
			{
				d00 = outside_value;
				d01 = outside_value;
			}
			else
			{
				d00 = A2D_ELEM(*this, y0, x0);
				if (outX1)
				{
					d01 = outside_value;
				}
				else
				{
					d01 = A2D_ELEM(*this, y0, x1);
				}
			}
		}
		else
		{
			if (outX0)
			{
				d00 = outside_value;
				d10 = outside_value;
				d01 = outside_value;
				d11 = outside_value;
			}
			else
			{
				if (outX1)
				{
					d01 = outside_value;
					d11 = outside_value;
				}
				else
				{
					d01 = A2D_ELEM(*this, y0, x1);
					d11 = A2D_ELEM(*this, y1, x1);
				}

				d00 = A2D_ELEM(*this, y0, x0);
				d10 = A2D_ELEM(*this, y1, x0);
			}
		}
    }

    double d0 = LIN_INTERP(fx, d00, d01);
    double d1 = LIN_INTERP(fx, d10, d11);
    return LIN_INTERP(fy, d0, d1);
}

void sincos(const MultidimArray<double> &x, MultidimArray<double> &s, MultidimArray<double> &c)
{
    s.resizeNoCopy(x);
    c.resizeNoCopy(x);
    double *ptr=NULL;
    double *ptrS=MULTIDIM_ARRAY(s);
    double *ptrC=MULTIDIM_ARRAY(c);
    size_t n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(x,n,ptr)
    sincos(*ptr, ptrS++,ptrC++);
}


void planeFit(const MultidimArray<double> &z, const MultidimArray<double> &x, const MultidimArray<double> &y,
		double &p0, double &p1, double &p2)
{
	 if (MULTIDIM_SIZE(z)!=MULTIDIM_SIZE(y) || MULTIDIM_SIZE(z)!=MULTIDIM_SIZE(x))
		 REPORT_ERROR(ERR_MULTIDIM_SIZE,"Not all vectors are of the same size");
	 if (MULTIDIM_SIZE(z) < 10)
		 REPORT_ERROR(ERR_MULTIDIM_SIZE, "Not enough elements to compute Least Squares plane fit");

	 double m11=0, m12=0, m13=0, m21=0, m22=0, m23=0, m31=0, m32=0, m33=0;
	 double b1=0, b2=0, b3=0;

	 FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(z)
	 {
		 double X=DIRECT_MULTIDIM_ELEM(x,n);
		 double Y=DIRECT_MULTIDIM_ELEM(y,n);
		 double Z=DIRECT_MULTIDIM_ELEM(z,n);
		 m11+=X*X;
		 m12+=X*Y;
		 m13+=X;

		 m22+=Y*Y;
		 m23+=Y;

		 b1+=X*Z;
		 b2+=Y*Z;
		 b3+=Z;
	 }
	 m21=m12;
	 m31=m13;
	 m32=m23;
	 m33=MULTIDIM_SIZE(z);

	 Matrix2D<double> A(3, 3);
	 Matrix1D<double> b(3);
	 Matrix1D<double> c(3);

	 A(0,0)=m11;
	 A(0,1)=m12;
	 A(0,2)=m13;
	 A(1,0)=m21;
	 A(1,1)=m22;
	 A(1,2)=m23;
	 A(2,0)=m31;
	 A(2,1)=m32;
	 A(2,2)=m33;

	 b(0)=b1;
	 b(1)=b2;
	 b(2)=b3;

	 c = A.inv() * b;
	 p0 = c(2);
	 p2 = c(1);
	 p1 = c(0);
}

// explicit instantiation
template MultidimArray<double>& MultidimArray<double>::operator=(Matrix2D<double> const&);
template void MultidimArray<double>::copy(Matrix2D<double>&) const;
template void MultidimArray<double>::rangeAdjustLS(const MultidimArray<double> &example, const MultidimArray<int> *mask);
