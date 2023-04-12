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

#include <fstream>
#include <algorithm>
#include "bilib/kernel.h"
#include "matrix2d.h"
#include "multidim_array.h"
#include "multidim_array_base.h"
#include "numerical_recipes.h"
#include "xmipp_funcs.h"
#include "xmipp_filename.h"
#include "utils/half.hpp"

template<typename T>
void MultidimArray<T>::getSliceAsMatrix(size_t k, Matrix2D<T> &m) const
{
//    there's wrong bracket somewhere, and I don't want to think about it now
//    m.resizeNoCopy(YSIZE(*this),XSIZE(*this));
//    memcpy(&MAT_ELEM(m,0,0),&A3D_ELEM(*this,k,0,0),YSIZE(*this),XSIZE(*this)*sizeof(double));
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Please contact developers");
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

template<typename T>
T MultidimArray<T>::interpolatedElementBSpline3D(double x, double y, double z,
                               int SplineDegree) const
{
    int SplineDegree_1 = SplineDegree - 1;

    // Logical to physical
    z -= STARTINGZ(*this);
    y -= STARTINGY(*this);
    x -= STARTINGX(*this);

    int l1 = (int)ceil(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = (int)ceil(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = (int)ceil(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    double aux;
    int Xdim=(int)XSIZE(*this);
    int Ydim=(int)YSIZE(*this);
    int Zdim=(int)ZSIZE(*this);
    for (int nn = n1; nn <= n2; nn++)
    {
        int equivalent_nn=nn;
        if      (nn<0)
            equivalent_nn=-nn-1;
        else if (nn>=Zdim)
            equivalent_nn=2*Zdim-nn-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=Ydim)
                equivalent_m=2*Ydim-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=Xdim)
                    equivalent_l=2*Xdim-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(*this,
                                                        equivalent_nn,equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    BSPLINE03(aux,xminusl);
                    xsum += Coeff * aux;
                    break;
                case 4:
                    xsum += Coeff * Bspline04(xminusl);
                    break;
                case 5:
                    xsum += Coeff * Bspline05(xminusl);
                    break;
                case 6:
                    xsum += Coeff * Bspline06(xminusl);
                    break;
                case 7:
                    xsum += Coeff * Bspline07(xminusl);
                    break;
                case 8:
                    xsum += Coeff * Bspline08(xminusl);
                    break;
                case 9:
                    xsum += Coeff * Bspline09(xminusl);
                    break;
                }
            }

            double yminusm = y - (double) m;
            switch (SplineDegree)
            {
            case 2:
                yxsum += xsum * Bspline02(yminusm);
                break;
            case 3:
                BSPLINE03(aux,yminusm);
                yxsum += xsum * aux;
                break;
            case 4:
                yxsum += xsum * Bspline04(yminusm);
                break;
            case 5:
                yxsum += xsum * Bspline05(yminusm);
                break;
            case 6:
                yxsum += xsum * Bspline06(yminusm);
                break;
            case 7:
                yxsum += xsum * Bspline07(yminusm);
                break;
            case 8:
                yxsum += xsum * Bspline08(yminusm);
                break;
            case 9:
                yxsum += xsum * Bspline09(yminusm);
                break;
            }
        }

        double zminusn = z - (double) nn;
        switch (SplineDegree)
        {
        case 2:
            zyxsum += yxsum * Bspline02(zminusn);
            break;
        case 3:
            BSPLINE03(aux,zminusn);
            zyxsum += yxsum * aux;
            break;
        case 4:
            zyxsum += yxsum * Bspline04(zminusn);
            break;
        case 5:
            zyxsum += yxsum * Bspline05(zminusn);
            break;
        case 6:
            zyxsum += yxsum * Bspline06(zminusn);
            break;
        case 7:
            zyxsum += yxsum * Bspline07(zminusn);
            break;
        case 8:
            zyxsum += yxsum * Bspline08(zminusn);
            break;
        case 9:
            zyxsum += yxsum * Bspline09(zminusn);
            break;
        }
    }

    return (T) zyxsum;
}

template<typename T>
T MultidimArray<T>::interpolatedElementBSpline2D(double x, double y, int SplineDegree) const
{
    int SplineDegree_1 = SplineDegree - 1;

    // Logical to physical
    y -= STARTINGY(*this);
    x -= STARTINGX(*this);

    int l1 = (int)ceil(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;
    int m1 = (int)ceil(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    double columns = 0.0;
    double aux;
    int Ydim=(int)YSIZE(*this);
    int Xdim=(int)XSIZE(*this);
    for (int m = m1; m <= m2; m++)
    {
        int equivalent_m=m;
        if      (m<0)
            equivalent_m=-m-1;
        else if (m>=Ydim)
            equivalent_m=2*Ydim-m-1;
        double rows = 0.0;
        for (int l = l1; l <= l2; l++)
        {
            double xminusl = x - (double) l;
            int equivalent_l=l;
            if      (l<0)
                equivalent_l=-l-1;
            else if (l>=Xdim)
                equivalent_l=2*Xdim-l-1;
            double Coeff = DIRECT_A2D_ELEM(*this, equivalent_m,equivalent_l);
            switch (SplineDegree)
            {
            case 2:
                rows += Coeff * Bspline02(xminusl);
                break;

            case 3:
                BSPLINE03(aux,xminusl);
                rows += Coeff * aux;
                break;

            case 4:
                rows += Coeff * Bspline04(xminusl);
                break;

            case 5:
                rows += Coeff * Bspline05(xminusl);
                break;

            case 6:
                rows += Coeff * Bspline06(xminusl);
                break;

            case 7:
                rows += Coeff * Bspline07(xminusl);
                break;

            case 8:
                rows += Coeff * Bspline08(xminusl);
                break;

            case 9:
                rows += Coeff * Bspline09(xminusl);
                break;
            }
        }

        double yminusm = y - (double) m;
        switch (SplineDegree)
        {
        case 2:
            columns += rows * Bspline02(yminusm);
            break;

        case 3:
            BSPLINE03(aux,yminusm);
            columns += rows * aux;
            break;

        case 4:
            columns += rows * Bspline04(yminusm);
            break;

        case 5:
            columns += rows * Bspline05(yminusm);
            break;

        case 6:
            columns += rows * Bspline06(yminusm);
            break;

        case 7:
            columns += rows * Bspline07(yminusm);
            break;

        case 8:
            columns += rows * Bspline08(yminusm);
            break;

        case 9:
            columns += rows * Bspline09(yminusm);
            break;
        }
    }
    return (T) columns;
}

template<typename T>
T MultidimArray<T>::interpolatedElementBSpline2D_Degree3(double x, double y) const
{
    bool    firstTime=true;         // Inner loop first time execution flag.
    double  *ref;

    // Logical to physical
    y -= STARTINGY(*this);
    x -= STARTINGX(*this);

    int l1 = (int)ceil(x - 2);
    int l2 = l1 + 3;
    int m1 = (int)ceil(y - 2);
    int m2 = m1 + 3;

    double columns = 0.0;
    double aux;
    int Ydim=(int)YSIZE(*this);
    int Xdim=(int)XSIZE(*this);

    int     equivalent_l_Array[LOOKUP_TABLE_LEN]; // = new int [l2 - l1 + 1];
    double  aux_Array[LOOKUP_TABLE_LEN];// = new double [l2 - l1 + 1];

    for (int m = m1; m <= m2; m++)
    {
        int equivalent_m=m;
        if      (m<0)
            equivalent_m=-m-1;
        else if (m>=Ydim)
            equivalent_m=2*Ydim-m-1;
        double rows = 0.0;
        int index=0;
        ref = &DIRECT_A2D_ELEM(*this, equivalent_m,0);
        for (int l = l1; l <= l2; l++)
        {
            int equivalent_l;
            // Check if it is first time executing inner loop.
            if (firstTime)
            {
                double xminusl = x - (double) l;
                equivalent_l=l;
                if (l<0)
                {
                    equivalent_l=-l-1;
                }
                else if (l>=Xdim)
                {
                    equivalent_l=2*Xdim-l-1;
                }

                equivalent_l_Array[index] = equivalent_l;
                BSPLINE03(aux,xminusl);
                aux_Array[index] = aux;
                index++;
            }
            else
            {
                equivalent_l = equivalent_l_Array[index];
                aux = aux_Array[index];
                index++;
            }

            //double Coeff = DIRECT_A2D_ELEM(*this, equivalent_m,equivalent_l);
            double Coeff = ref[equivalent_l];
            rows += Coeff * aux;
        }

        // Set first time inner flag is executed to false.
        firstTime = false;

        double yminusm = y - (double) m;
        BSPLINE03(aux,yminusm);
        columns += rows * aux;
    }

    return (T) columns;
}

template<typename T>
T MultidimArray<T>::interpolatedElementBSpline1D(double x, int SplineDegree) const
{
    int SplineDegree_1 = SplineDegree - 1;

    // Logical to physical
    x -= STARTINGX(*this);

    int l1 = (int)ceil(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;
    int Xdim=(int)XSIZE(*this);
    double sum = 0.0;
    for (int l = l1; l <= l2; l++)
    {
        double xminusl = x - (double) l;
        int equivalent_l=l;
        if      (l<0)
            equivalent_l=-l-1;
        else if (l>=Xdim)
            equivalent_l=2*Xdim-l-1;
        double Coeff = (double) DIRECT_A1D_ELEM(*this, equivalent_l);
        double aux;
        switch (SplineDegree)
        {
        case 2:
            sum += Coeff * Bspline02(xminusl);
            break;

        case 3:
            BSPLINE03(aux,xminusl);
            sum += Coeff * aux;
            break;

        case 4:
            sum += Coeff * Bspline04(xminusl);
            break;

        case 5:
            sum += Coeff * Bspline05(xminusl);
            break;

        case 6:
            sum += Coeff * Bspline06(xminusl);
            break;

        case 7:
            sum += Coeff * Bspline07(xminusl);
            break;

        case 8:
            sum += Coeff * Bspline08(xminusl);
            break;

        case 9:
            sum += Coeff * Bspline09(xminusl);
            break;
        }
    }
    return (T) sum;
}

template<typename T>
void MultidimArray<T>::initRandom(double op1, double op2, RandomMode mode)
{
    T* ptr=NULL;
    size_t n;
    if (mode == RND_UNIFORM)
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = static_cast< T >(rnd_unif(op1, op2));
    else if (mode == RND_GAUSSIAN)
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = static_cast< T >(rnd_gaus(op1, op2));
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT,
                     formatString("InitRandom: Mode not supported"));
}

template<typename T>
void MultidimArray<T>::addNoise(double op1,
              double op2,
              const String& mode,
              double df) const
{
    T* ptr=NULL;
    size_t n;
    if (mode == "uniform")
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr += static_cast< T >(rnd_unif(op1, op2));
    else if (mode == "gaussian")
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr += static_cast< T >(rnd_gaus(op1, op2));
    else if (mode == "student")
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr += static_cast< T >(rnd_student_t(df, op1, op2));
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT,
                     formatString("AddNoise: Mode not supported (%s)", mode.c_str()));
}

template<typename T>
FILE* MultidimArray<T>::mmapFile(T* &_data, size_t nzyxDim) const
{
#ifdef XMIPP_MMAP
    FILE* fMap = tmpfile();
    int Fd = fileno(fMap);

    if ((lseek(Fd, nzyxDim*sizeof(T)-1, SEEK_SET) == -1) || (::write(Fd,"",1) == -1))
    {
        fclose(fMap);
        REPORT_ERROR(ERR_IO_NOWRITE,"MultidimArray::resize: Error 'stretching' the map file.");
    }
    if ( (_data = (T*) mmap(0,nzyxDim*sizeof(T), PROT_READ | PROT_WRITE, MAP_SHARED, Fd, 0)) == (void*) MAP_FAILED )
        REPORT_ERROR(ERR_MMAP_NOTADDR,formatString("MultidimArray::resize: mmap failed. Error %s", strerror(errno)));

    return fMap;
#else

    REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

}

template<typename T>
void MultidimArray<T>::indexSort(MultidimArray< int > &indx) const
{
    checkDimension(1);

    MultidimArray< double > temp;
    indx.clear();

    if (xdim == 0)
        return;

    if (xdim == 1)
    {
        indx.resizeNoCopy(1);
        DIRECT_A1D_ELEM(indx,0) = 1;
        return;
    }

    // Initialise data
    indx.resizeNoCopy(xdim);
    typeCast(*this, temp);

    // Sort indexes
    double* temp_array = temp.adaptForNumericalRecipes1D();
    int* indx_array = indx.adaptForNumericalRecipes1D();
    indexx(XSIZE(*this), temp_array, indx_array);
}

template<typename T>
void MultidimArray<T>::selfNormalizeInterval(double minPerc, double maxPerc, int Npix)
{
    std::vector<double> randValues; // Vector with random chosen values

    for(int i=0; i<Npix; i++)
    {
        size_t indx = (size_t)rnd_unif(0, MULTIDIM_SIZE(*this));
        randValues.push_back(DIRECT_MULTIDIM_ELEM(*this,indx));
    }
    std::sort(randValues.begin(),randValues.end());

    double m = randValues[(size_t)(minPerc*Npix)];
    double M = randValues[(size_t)(maxPerc*Npix)];

    T* ptr=NULL;
    size_t n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
    *ptr = 2/(M-m)*(*ptr-m)-1;
}

template<typename T>
void MultidimArray<T>::showWithGnuPlot(const String& xlabel, const String& title)
{
    checkDimension(1);

    FileName fn_tmp;
    fn_tmp.initRandom(10);
    const char * fnStr = fn_tmp.c_str();
    MultidimArray<T>::write(formatString("PPP%s.txt", fnStr));

    std::ofstream fh_gplot;
    fh_gplot.open(formatString("PPP%s.gpl", fnStr).c_str());
    if (!fh_gplot)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     formatString("vector::showWithGnuPlot: Cannot open PPP%s.gpl for output", fnStr));
    fh_gplot << "set xlabel \"" + xlabel + "\"\n";
    fh_gplot << "plot \"PPP" + fn_tmp + ".txt\" title \"" + title +
    "\" w l\n";
    fh_gplot << "pause 300 \"\"\n";
    fh_gplot.close();
    if (0 != system(formatString("(gnuplot PPP%s.gpl; rm PPP%s.txt PPP%s.gpl) &", fnStr, fnStr, fnStr).c_str()) ) {
        REPORT_ERROR(ERR_IO, "Cannot open gnuplot");
    }
}

template<typename T>
void MultidimArray<T>::edit()
{
    FileName nam;
    nam.initRandom(15);

    nam = formatString("PPP%s.txt", nam.c_str());
    write(nam);

    if (0 != system(formatString("xmipp_edit -i %s -remove &", nam.c_str()).c_str())) {
        REPORT_ERROR(ERR_IO, "Cannot open xmipp_edit");
    }
}

template<typename T>
void MultidimArray<T>::write(const FileName& fn) const
{
    std::ofstream out;
    out.open(fn.c_str(), std::ios::out);
    if (!out)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     formatString("MultidimArray::write: File %s cannot be opened for output", fn.c_str()));

    out << *this;
    out.close();
}

template<typename T>
void MultidimArray<T>::resize(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim, bool copy)
{
    if (Ndim*Zdim*Ydim*Xdim == nzyxdimAlloc && data != NULL)
    {
        ndim = Ndim;
        xdim = Xdim;
        ydim = Ydim;
        zdim = Zdim;
        yxdim = Ydim * Xdim;
        zyxdim = Zdim * yxdim;
        nzyxdim = Ndim * zyxdim;
        return;
    }
    else if (!destroyData)
        REPORT_ERROR(ERR_MULTIDIM_SIZE, "Cannot resize array when accessing through alias.");

    if (Xdim <= 0 || Ydim <= 0 || Zdim <= 0 || Ndim <= 0)
    {
        clear();
        return;
    }

    // data can be NULL while xdim etc are set to non-zero values
    // (This can happen for reading of images...)
    // In that case, initialize data to zeros.
    if (NZYXSIZE(*this) > 0 && data == NULL)
    {
        ndim = Ndim;
        xdim = Xdim;
        ydim = Ydim;
        zdim = Zdim;
        yxdim = Ydim * Xdim;
        zyxdim = Zdim * yxdim;
        nzyxdim = Ndim * zyxdim;

        coreAllocate();
        return;
    }

    // Ask for memory
    size_t YXdim=(size_t)Ydim*Xdim;
    size_t ZYXdim=YXdim*Zdim;
    size_t NZYXdim=ZYXdim*Ndim;
    FILE*  new_mFd=NULL;

    T * new_data=NULL;

    try
    {
        if (mmapOn)
            new_mFd = mmapFile(new_data, NZYXdim);
        else
            new_data = new T [NZYXdim];

        memset(new_data,0,NZYXdim*sizeof(T));
    }
    catch (std::bad_alloc &)
    {
        if (!mmapOn)
        {
            setMmap(true);
            resize(Ndim, Zdim, Ydim, Xdim, copy);
            return;
        }
        else
        {
            std::ostringstream sstream;
            sstream << "Allocate: No space left to allocate ";
            sstream << (NZYXdim * sizeof(T)/1024/1024/1024) ;
            sstream << "Gb." ;
            REPORT_ERROR(ERR_MEM_NOTENOUGH, sstream.str());
        }
    }
    // Copy needed elements, fill with 0 if necessary
    if (copy)
    {
        const auto nCopy = std::min(Ndim, NSIZE(*this));
        const auto zCopy = std::min(Zdim, ZSIZE(*this));
        const auto yCopy = std::min(Ydim, YSIZE(*this));
        const auto xCopy = std::min(Xdim, XSIZE(*this));
        const auto xCopyBytes = xCopy*sizeof(T);

        for (size_t l = 0; l < nCopy; ++l) 
            for (size_t k = 0; k < zCopy; ++k) 
                for (size_t i = 0; i < yCopy; ++i) 
                    memcpy(
                        new_data + l*ZYXdim + k*YXdim + i*Xdim,
                        data + l*zyxdim + k*yxdim + i*xdim,
                        xCopyBytes
                    );
    }

    // deallocate old array
    coreDeallocate();

    // assign *this vector to the newly created
    data = new_data;
    ndim = Ndim;
    xdim = Xdim;
    ydim = Ydim;
    zdim = Zdim;
    yxdim = Ydim * Xdim;
    zyxdim = Zdim * yxdim;
    nzyxdim = Ndim * zyxdim;
    mFd = new_mFd;
    nzyxdimAlloc = nzyxdim;
}

template<typename T>
void MultidimArray<T>::sort(MultidimArray<T> &result) const
{
    result = *this;
    std::sort(result.data, result.data + result.nzyxdim);
}

template<typename T>
void MultidimArray<T>::computeMedian_within_binary_mask(const MultidimArray< int >& mask, double& median) const
{
    std::vector<double> bgI;

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*this)
    {
        if (DIRECT_MULTIDIM_ELEM(mask, n) != 0)
        {
            double aux = DIRECT_MULTIDIM_ELEM(*this, n);
            bgI.push_back(aux);
        }
    }

    std::sort(bgI.begin(), bgI.end());
    if (bgI.size() % 2 != 0)
        median = bgI[bgI.size() / 2];
    else
        median = (bgI[(bgI.size() - 1) / 2] + bgI[bgI.size() / 2]) / 2.0;
}

template<typename T>
void MultidimArray<T>::randomSubstitute(T oldv,
                      T avgv,
                      T sigv,
                      double accuracy,
                      MultidimArray<int> * mask)
{
    T* ptr=NULL;
    size_t n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
    if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
        if (ABS(*ptr - oldv) <= accuracy)
            *ptr = rnd_gaus(avgv, sigv);
}

// explicit instantiation
template class MultidimArray<double>;

// mmapFile
template FILE* MultidimArray<bool>::mmapFile(bool*&, unsigned long) const;
template FILE* MultidimArray<float>::mmapFile(float*&, unsigned long) const;
template FILE* MultidimArray<char>::mmapFile(char*&, unsigned long) const;
template FILE* MultidimArray<int>::mmapFile(int*&, unsigned long) const;
template FILE* MultidimArray<long>::mmapFile(long*&, unsigned long) const;
template FILE* MultidimArray<short>::mmapFile(short*&, unsigned long) const;
template FILE* MultidimArray<unsigned char>::mmapFile(unsigned char*&, unsigned long) const;
template FILE* MultidimArray<unsigned int>::mmapFile(unsigned int*&, unsigned long) const;
template FILE* MultidimArray<unsigned long>::mmapFile(unsigned long*&, unsigned long) const;
template FILE* MultidimArray<unsigned short>::mmapFile(unsigned short*&, unsigned long) const;
template FILE* MultidimArray<std::complex<double>>::mmapFile(std::complex<double>*&, unsigned long) const;
template FILE* MultidimArray<half_float::half>::mmapFile(half_float::half*&, unsigned long) const;

// resize
template void MultidimArray<bool>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<float>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<char>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<long>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<int>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<short>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<unsigned char>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<unsigned int>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<unsigned long>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<unsigned short>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<std::complex<double>>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);
template void MultidimArray<half_float::half>::resize(unsigned long, unsigned long, unsigned long, unsigned long, bool);

// index sort
template void MultidimArray<float>::indexSort(MultidimArray<int>&) const;
template void MultidimArray<int>::indexSort(MultidimArray<int>&) const;

// init random
template void MultidimArray<float>::initRandom(double, double, RandomMode);
template void MultidimArray<char>::initRandom(double, double, RandomMode);
template void MultidimArray<long>::initRandom(double, double, RandomMode);
template void MultidimArray<int>::initRandom(double, double, RandomMode);
template void MultidimArray<short>::initRandom(double, double, RandomMode);
template void MultidimArray<unsigned char>::initRandom(double, double, RandomMode);
template void MultidimArray<unsigned int>::initRandom(double, double, RandomMode);
template void MultidimArray<unsigned long>::initRandom(double, double, RandomMode);
template void MultidimArray<unsigned short>::initRandom(double, double, RandomMode);
template void MultidimArray<half_float::half>::initRandom(double, double, RandomMode);

// write
template void MultidimArray<bool>::write(FileName const&) const;
