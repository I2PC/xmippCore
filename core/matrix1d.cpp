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

#include "matrix1d.h"
#include "xmipp_filename.h"
#include "numerical_recipes.h"
#include <fstream>
#include <algorithm>

template<typename T>
Matrix1D<T>& Matrix1D<T>::operator=(const Matrix1D<T>& op1)
{
    if (&op1 != this)
    {
        resizeNoCopy(op1);
        memcpy(vdata, op1.vdata, vdim * sizeof(T));
        row = op1.row;
    }

    return *this;
}

template<typename T>
bool Matrix1D<T>::operator==(const Matrix1D<T>& op1) const
{
    if (row != op1.row)
        return false;

    for (size_t i = 0; i < vdim; ++i)
        if (!XMIPP_EQUAL_REAL(vdata[i], op1.vdata[i]))
            return false;
    return true;
}

template<typename T>
Matrix1D<T>& Matrix1D<T>::operator=(const std::vector<T>& op1)
{
    resizeNoCopy(op1.size());
    memcpy(&vdata[0],&(op1[0]),vdim*sizeof(T));
    row = false;
    return *this;
}

template<typename T>
void Matrix1D<T>::clear()
{
    coreDeallocate();
    coreInit();
}

template<typename T>
void Matrix1D<T>::coreInit()
{
    vdim = 0;
    row = false;
    vdata = NULL;
    destroyData = true;
}

template<typename T>
void Matrix1D<T>::initConstant(T val)
{
    for (size_t j = 0; j < vdim; j++)
        vdata[j] = val;
}

template<typename T>
void Matrix1D<T>::initConstant(size_t Xdim, T val)
{
    if (vdim != Xdim)
        resizeNoCopy(Xdim);
    initConstant(val);
}

template<typename T>
void Matrix1D<T>::enumerate()
{
  for (size_t j = 0; j < vdim; j++)
      vdata[j] = j;
}

template<typename T>
bool Matrix1D<T>::isAnyNaN()
{
    for (size_t j = 0; j < vdim; j++)
        if (ISNAN(vdata[j]))
            return true;
    return false;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(T op1) const
{
    Matrix1D<T> tmp(*this);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) * op1;
        (*ptr1++) = (*ptr2++) * op1;
        (*ptr1++) = (*ptr2++) * op1;
        (*ptr1++) = (*ptr2++) * op1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) * op1;
    return tmp;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator/(T op1) const
{
    Matrix1D<T> tmp(*this);
    T iop1 = 1 / op1;
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) * iop1;
        (*ptr1++) = (*ptr2++) * iop1;
        (*ptr1++) = (*ptr2++) * iop1;
        (*ptr1++) = (*ptr2++) * iop1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) * iop1;
    return tmp;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator+(T op1) const
{
    Matrix1D<T> tmp(*this);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) + op1;
        (*ptr1++) = (*ptr2++) + op1;
        (*ptr1++) = (*ptr2++) + op1;
        (*ptr1++) = (*ptr2++) + op1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) + op1;
    return tmp;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator-(T op1) const
    {
        Matrix1D<T> tmp(*this);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) - op1;
            (*ptr1++) = (*ptr2++) - op1;
            (*ptr1++) = (*ptr2++) - op1;
            (*ptr1++) = (*ptr2++) - op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) - op1;
        return tmp;
    }

template<typename T>
void Matrix1D<T>::operator+=(const Matrix1D<T>& op1) const
{
    if (vdim != op1.vdim)
        REPORT_ERROR(ERR_MATRIX_SIZE, "Not same sizes in vector summation");

    T *ptr1 = &VEC_ELEM(*this,0);
    const T *ptr2 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) += (*ptr2++);
        (*ptr1++) += (*ptr2++);
        (*ptr1++) += (*ptr2++);
        (*ptr1++) += (*ptr2++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) += (*ptr2++);
}

template<typename T>
void Matrix1D<T>::operator-=(const Matrix1D<T>& op1) const
{
    if (vdim != op1.vdim)
        REPORT_ERROR(ERR_MATRIX_SIZE, "Not same sizes in vector summation");

    T *ptr1 = &VEC_ELEM(*this,0);
    const T *ptr2 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) -= (*ptr2++);
        (*ptr1++) -= (*ptr2++);
        (*ptr1++) -= (*ptr2++);
        (*ptr1++) -= (*ptr2++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) -= (*ptr2++);
}

template<typename T>
void Matrix1D<T>::operator*=(T op1)
{
    T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) *= op1;
        (*ptr1++) *= op1;
        (*ptr1++) *= op1;
        (*ptr1++) *= op1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) *= op1;
}

template<typename T>
void Matrix1D<T>::operator/=(T op1)
{
    T iop1 = 1 / op1;
    T * ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) *= iop1;
        (*ptr1++) *= iop1;
        (*ptr1++) *= iop1;
        (*ptr1++) *= iop1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) *= iop1;
}

template<typename T>
void Matrix1D<T>::operator+=(T op1)
{
    T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) += op1;
        (*ptr1++) += op1;
        (*ptr1++) += op1;
        (*ptr1++) += op1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) += op1;
}

template<typename T>
void Matrix1D<T>::operator-=(T op1)
{
    T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) -= op1;
        (*ptr1++) -= op1;
        (*ptr1++) -= op1;
        (*ptr1++) -= op1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) -= op1;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(const Matrix1D<T>& op1) const
{
    Matrix1D<T> tmp(op1);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    const T *ptr3 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) * (*ptr3++);
        (*ptr1++) = (*ptr2++) * (*ptr3++);
        (*ptr1++) = (*ptr2++) * (*ptr3++);
        (*ptr1++) = (*ptr2++) * (*ptr3++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) * (*ptr3++);
    return tmp;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator/(const Matrix1D<T>& op1) const
{
    Matrix1D<T> tmp(op1);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    const T *ptr3 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) / (*ptr3++);
        (*ptr1++) = (*ptr2++) / (*ptr3++);
        (*ptr1++) = (*ptr2++) / (*ptr3++);
        (*ptr1++) = (*ptr2++) / (*ptr3++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) / (*ptr3++);
    return tmp;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator+(const Matrix1D<T>& op1) const
{
    Matrix1D<T> tmp(op1);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    const T *ptr3 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) + (*ptr3++);
        (*ptr1++) = (*ptr2++) + (*ptr3++);
        (*ptr1++) = (*ptr2++) + (*ptr3++);
        (*ptr1++) = (*ptr2++) + (*ptr3++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) + (*ptr3++);
    return tmp;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator-(const Matrix1D<T>& op1) const
{
    Matrix1D<T> tmp(op1);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    const T *ptr3 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) = (*ptr2++) - (*ptr3++);
        (*ptr1++) = (*ptr2++) - (*ptr3++);
        (*ptr1++) = (*ptr2++) - (*ptr3++);
        (*ptr1++) = (*ptr2++) - (*ptr3++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) = (*ptr2++) - (*ptr3++);
    return tmp;
}

template<typename T>
void Matrix1D<T>::operator*=(const Matrix1D<T>& op1)
{
    T *ptr1 = &VEC_ELEM(*this,0);
    const T *ptr2 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) *= (*ptr2++);
        (*ptr1++) *= (*ptr2++);
        (*ptr1++) *= (*ptr2++);
        (*ptr1++) *= (*ptr2++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) *= (*ptr2++);
}

template<typename T>
void Matrix1D<T>::operator/=(const Matrix1D<T>& op1)
{
    T *ptr1 = &VEC_ELEM(*this,0);
    const T *ptr2 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        (*ptr1++) /= (*ptr2++);
        (*ptr1++) /= (*ptr2++);
        (*ptr1++) /= (*ptr2++);
        (*ptr1++) /= (*ptr2++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        (*ptr1++) /= (*ptr2++);
}

template<typename T>
void Matrix1D<T>::operator+=(const Matrix1D<T>& op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) += (*ptr2++);
    }

template<typename T>
void Matrix1D<T>::operator-=(const Matrix1D<T>& op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) -= (*ptr2++);
    }

template<typename T>
Matrix1D<T> Matrix1D<T>::operator-() const
{
    Matrix1D<T> tmp(*this);
    T *ptr1 = &VEC_ELEM(tmp,0);
    const T *ptr2 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        if constexpr (std::is_signed_v<T> || std::is_floating_point_v<T>) {
            (*ptr1++) = -(*ptr2++);
            (*ptr1++) = -(*ptr2++);
            (*ptr1++) = -(*ptr2++);
            (*ptr1++) = -(*ptr2++);
        } else {
            REPORT_ERROR(
            ERR_TYPE_INCORRECT,
            static_cast<std::string>("Can not use unitary minus on unsigned datatypes"));
        }
        
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        if constexpr (std::is_signed_v<T> || std::is_floating_point_v<T>) {
            (*ptr1++) = -(*ptr2++);
        } else {
            REPORT_ERROR(
            ERR_TYPE_INCORRECT,
            static_cast<std::string>("Can not use unitary minus on unsigned datatypes"));
        }
    return tmp;
}

template<typename T>
void Matrix1D<T>::selfCEIL()
{
    T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        *ptr1 = ceil(*ptr1);
        ++ptr1;
        *ptr1 = ceil(*ptr1);
        ++ptr1;
        *ptr1 = ceil(*ptr1);
        ++ptr1;
        *ptr1 = ceil(*ptr1);
        ++ptr1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
    {
        *ptr1 = ceil(*ptr1);
        ++ptr1;
    }
}

template<typename T>
void Matrix1D<T>::selfFLOOR()
{
    T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        *ptr1 = floor(*ptr1);
        ++ptr1;
        *ptr1 = floor(*ptr1);
        ++ptr1;
        *ptr1 = floor(*ptr1);
        ++ptr1;
        *ptr1 = floor(*ptr1);
        ++ptr1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
    {
        *ptr1 = floor(*ptr1);
        ++ptr1;
    }
}

template<typename T>
void Matrix1D<T>::selfROUND()
{
    T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        *ptr1 = round(*ptr1);
        ++ptr1;
        *ptr1 = round(*ptr1);
        ++ptr1;
        *ptr1 = round(*ptr1);
        ++ptr1;
        *ptr1 = round(*ptr1);
        ++ptr1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
    {
        *ptr1 = round(*ptr1);
        ++ptr1;
    }
}

    /** Mean value of the vector */
template<typename T>
double Matrix1D<T>::computeMean() const
{
    if (vdim == 0)
        return 0;

    double sum = 0;
    for (size_t j = 0; j < vdim; ++j)
        sum+=VEC_ELEM(*this,j);
    return sum/vdim;
}

template<typename T>
T Matrix1D<T>::computeMax() const
{
    if (vdim == 0)
        return 0;

    T maxval = VEC_ELEM(*this,0);
    for (size_t j = 0; j < vdim; ++j)
        if (VEC_ELEM(*this,j) > maxval)
            maxval = VEC_ELEM(*this,j);
    return maxval;
}

template<typename T>
void Matrix1D<T>::computeMinMax(T &minval, T &maxval) const
{
    if (vdim == 0)
        return;

    maxval = minval = VEC_ELEM(*this,0);
    for (size_t j = 0; j < vdim; ++j)
    {
        T val=VEC_ELEM(*this,j);
        if (val > maxval)
            maxval = val;
        else if (val<minval)
            minval = val;
    }
}

template<typename T>
void Matrix1D<T>::computeMeanAndStddev(double &mean, double &stddev) const
{
    mean=stddev=0;
    if (vdim == 0)
        return;

    double sum = 0, sum2 = 0;
    for (size_t j = 0; j < vdim; ++j)
    {
        double val=VEC_ELEM(*this,j);
        sum+=val;
        sum2+=val*val;
    }
    mean=sum/vdim;
    stddev=sum2/vdim-mean*mean;
    if (stddev<0)
        stddev=0;
    else
        stddev=sqrt(stddev);
}

template<typename T>
int Matrix1D<T>::maxIndex() const
{
    if (vdim == 0)
        return -1;

    int jmax = 0;
    T maxval = VEC_ELEM(*this, 0);
    for (size_t j = 0; j < vdim; ++j)
        if (VEC_ELEM(*this,j) > maxval)
        {
            jmax = j;
            maxval = VEC_ELEM(*this,j);
        }
    return jmax;
}

template<typename T>
int Matrix1D<T>::minIndex() const
{
    if (vdim == 0)
        return -1;

    int jmin = 0;
    T minval = VEC_ELEM(*this, 0);
    for (size_t j = 0; j < vdim; ++j)
        if (VEC_ELEM(*this,j) < minval)
        {
            jmin = j;
            minval = VEC_ELEM(*this,j);
        }
    return jmin;
}

template<typename T>
Matrix1D<T> Matrix1D<T>::transpose() const
{
    Matrix1D<T> temp(*this);
    temp.selfTranspose();
    return temp;
}

template<typename T>
double Matrix1D<T>::sum(bool average) const
{
    double sum = 0;
    const T *ptr1 = &VEC_ELEM(*this,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        sum += (*ptr1++);
        sum += (*ptr1++);
        sum += (*ptr1++);
        sum += (*ptr1++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        sum += (*ptr1++);
    if (average)
        return sum / (double) vdim;
    else
        return sum;
}

template<typename T>
double Matrix1D<T>::sum2() const
{
    double sum = 0;
    const T *ptr1 = &VEC_ELEM(*this,0);
    double val;
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        val = *ptr1;
        sum += val * val;
        ++ptr1;
        val = *ptr1;
        sum += val * val;
        ++ptr1;
        val = *ptr1;
        sum += val * val;
        ++ptr1;
        val = *ptr1;
        sum += val * val;
        ++ptr1;
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
    {
        val = *ptr1;
        sum += val * val;
        ++ptr1;
    }
    return sum;
}

template<typename T>
double Matrix1D<T>::dotProduct(const Matrix1D<T> &op1) const
{
    double sum = 0;
    const T *ptr1 = &VEC_ELEM(*this,0);
    const T *ptr2 = &VEC_ELEM(op1,0);
    size_t iBlockMax = vdim / 4;
    for (size_t i = 0; i < iBlockMax; i++)
    {
        sum += (*ptr1++) * (*ptr2++);
        sum += (*ptr1++) * (*ptr2++);
        sum += (*ptr1++) * (*ptr2++);
        sum += (*ptr1++) * (*ptr2++);
    }
    for (size_t i = iBlockMax * 4; i < vdim; ++i)
        sum += (*ptr1++) * (*ptr2++);
    return sum;
}

template<typename T>
void Matrix1D<T>::selfNormalize()
{
    double m = module();
    if (fabs(m) > XMIPP_EQUAL_ACCURACY)
    {
        T im = (T) (1.0 / m);
        *this *= im;
    }
    else
        initZeros();
}

template<typename T>
void Matrix1D<T>::selfReverse()
{
    for (int j = 0; j <= (int) (vdim - 1) / 2; j++)
    {
        T aux;
        SWAP(vdata[j], vdata[vdim-1-j], aux);
    }
}

template<typename T>
void Matrix1D<T>::numericalDerivative(Matrix1D<double> &result) const
{
    const double i12 = 1.0 / 12.0;
    result.initZeros(*this);
    for (int i = 2; i <= vdim - 2; i++)
        if constexpr (std::is_signed_v<T> || std::is_floating_point_v<T>) {
            result(i) = i12
                    * (-(*this)(i + 2) + 8 * (*this)(i + 1) - 8 * (*this)(i - 1)
                       + (*this)(i + 2));
        } else {
            REPORT_ERROR(
            ERR_TYPE_INCORRECT,
            static_cast<std::string>("Can not use unitary minus on unsigned datatypes"));
        }
}

template<typename T>
void Matrix1D<T>::showWithGnuPlot(const std::string& xlabel, const std::string& title)
{
    FileName fn_tmp;
    fn_tmp.initRandom(10);
    Matrix1D<T>::write(static_cast<std::string>("PPP") + fn_tmp + ".txt");

    std::ofstream fh_gplot;
    fh_gplot.open(
        (static_cast<std::string>("PPP") + fn_tmp + ".gpl").c_str());
    if (!fh_gplot)
        REPORT_ERROR(
            ERR_UNCLASSIFIED,
            static_cast<std::string>("vector::showWithGnuPlot: Cannot open PPP") + fn_tmp + ".gpl for output");
    fh_gplot << "set xlabel \"" + xlabel + "\"\n";
    fh_gplot
    << "plot \"PPP" + fn_tmp + ".txt\" title \"" + title
    + "\" w l\n";
    fh_gplot << "pause 300 \"\"\n";
    fh_gplot.close();
    auto res = system(
        (static_cast<std::string>("(gnuplot PPP") + fn_tmp
         + ".gpl; rm PPP" + fn_tmp + ".txt PPP" + fn_tmp
         + ".gpl) &").c_str());
    if (0 != res) {
        REPORT_ERROR(
            ERR_UNCLASSIFIED,
            "Something went wrong when working with GNUPlot. Please report this error to developers.");
    }
}

template<typename T>
void Matrix1D<T>::write(const FileName& fn) const
{
    std::ofstream out;
    out.open(fn.c_str(), std::ios::out);
    if (!out)
        REPORT_ERROR(
            ERR_IO_NOTOPEN,
            static_cast< std::string >("Matrix1D::write: File " + fn + " cannot be opened for output"));

    out << *this;
    out.close();
}

template<typename T>
void Matrix1D<T>::read(const FileName& fn)
{
    std::ifstream in;
    in.open(fn.c_str(), std::ios::in);

    if (!in)
        REPORT_ERROR(
            ERR_IO_NOTOPEN,
            static_cast< std::string >("MultidimArray::read: File " + fn + " not found"));

    in >> *this;
    in.close();
}

template<typename T>
void Matrix1D<T>::edit()
{
    FileName nam;
    nam.initRandom(15);

    nam = static_cast<std::string>("PPP" + nam + ".txt");
    write
    (nam);

    auto res = system(
        (static_cast<std::string>("xmipp_edit -i " + nam + " -remove &").c_str()));
    if (0 != res) {
        REPORT_ERROR(
            ERR_UNCLASSIFIED,
            "Something went wrong when working with xmipp_edit. Please report this error to developers.");
    }
}

template<typename T>
T Matrix1D<T>::computeMedian() const
{
    Matrix1D<T> aux;
    aux=*this;
    std::nth_element(aux.vdata,aux.vdata+vdim/2,aux.vdata+vdim);
    return VEC_ELEM(aux,vdim/2);
}
template<typename T>
Matrix1D<T> Matrix1D<T>::sort() const
{
    Matrix1D<T> temp(*this);
    if (vdim == 0)
        return temp;

    std::sort(temp.vdata, temp.vdata + temp.vdim);
    return temp;
}

template<typename T>
void Matrix1D<T>::indexSort(Matrix1D<int> &indx) const
{
    Matrix1D< double > temp;
    indx.clear();

    if (VEC_XSIZE(*this) == 0)
        return;

    if (VEC_XSIZE(*this) == 1)
    {
        indx.resizeNoCopy(1);
        VEC_ELEM(indx,0) = 1;
        return;
    }

    // Initialise data
    indx.resizeNoCopy(VEC_XSIZE(*this));
    typeCast(*this, temp);

    // Sort indexes
    indexx(VEC_XSIZE(*this), MATRIX1D_ARRAY(temp)-1, MATRIX1D_ARRAY(indx)-1);
}

Matrix1D<double> vectorR2(double x, double y)
{
    Matrix1D<double> result(2);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    return result;
}

Matrix1D<double> vectorR3(double x, double y, double z)
{
    Matrix1D<double> result(3);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    VEC_ELEM(result, 2) = z;
    return result;
}

Matrix1D<int> vectorR3(int x, int y, int z)
{
    Matrix1D<int> result(3);
    VEC_ELEM(result, 0) = x;
    VEC_ELEM(result, 1) = y;
    VEC_ELEM(result, 2) = z;
    return result;
}

// explicit instantiation
template class Matrix1D<double>;
template class Matrix1D<int>;
template class Matrix1D<float>;
template class Matrix1D<unsigned char>;
template class Matrix1D<short>;
template class Matrix1D<unsigned long>;
