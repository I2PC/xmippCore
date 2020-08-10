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
#include <fstream>

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
    system(
        (static_cast<std::string>("(gnuplot PPP") + fn_tmp
         + ".gpl; rm PPP" + fn_tmp + ".txt PPP" + fn_tmp
         + ".gpl) &").c_str());
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

    system(
        (static_cast<std::string>("xmipp_edit -i " + nam + " -remove &").c_str()));
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
