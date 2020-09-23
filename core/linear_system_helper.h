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

#ifndef XMIPPCORE_CORE_LINEAR_SYSTEM_HELPER_H_
#define XMIPPCORE_CORE_LINEAR_SYSTEM_HELPER_H_

#include "matrix1d.h"
#include "matrix2d.h"

/** Helper class for solving linear systems */
class PseudoInverseHelper
{
public:
    Matrix2D<double> A, AtA, AtAinv;
    Matrix1D<double> Atb, b, bpredicted;
};

/** Helper class for solving Weighted Least Squares */
class WeightedLeastSquaresHelper: public PseudoInverseHelper
{
public:
    Matrix1D<double> w; //Weights
};

/** Solve Linear system Ax=b with pseudoinverse.
 * A and b must be set inside the PseudoInverseHelper, the rest of the
 * fields in PseudoInverseHelper are used by this routine to avoid
 * several allocation/deallocations */
void solveLinearSystem(PseudoInverseHelper &h, Matrix1D<double> &result);

/** Solve Weighted least square problem Ax=b with pseudoinverse and weights w.
 * A, w and b must be set inside the WeightedLeastSquaresHelper, the rest of the
 * fields in WeightedLeastSquaresHelper are used by this routine to avoid
 * several allocation/deallocations.
 *
 * The normal equations of this problem are A^t W A x = A^t W b,
 * where W is a diagonal matrix whose entries are in the vector w. */
void weightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result);

/** Solve Weighted least square problem Ax=b and weights w, with RANSAC.
 * Tol is a tolerance value: if an equation is fulfilled with an error smaller than tol,
 * then it is considered to fulfill the model. Niter is the number of RANSAC iterations to perform.
 * The outlier fraction is the fraction of equations that, at maximum, can be considered as outliers.
 */
void ransacWeightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result,
        double tol, int Niter=10000, double outlierFraction=0.25, int Nthreads=1);

#endif /* XMIPPCORE_CORE_LINEAR_SYSTEM_HELPER_H_ */