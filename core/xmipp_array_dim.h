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

#ifndef XMIPPCORE_CORE_XMIPP_ARRAY_DIM_H_
#define XMIPPCORE_CORE_XMIPP_ARRAY_DIM_H_

#include <stddef.h>

/**
 *  Structure with the dimensions information of an image
 */
struct ArrayDim
{
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

    ArrayDim()
    {
        ndim = 0;
        zdim = 0;
        ydim = 0;
        xdim = 0;
        yxdim = 0;
        zyxdim = 0;
        nzyxdim = 0;
    }

    bool operator==(ArrayDim &adim)
    {
        return (this->ndim == adim.ndim &&
                this->zdim == adim.zdim &&
                this->ydim == adim.ydim &&
                this->xdim == adim.xdim );
    }
}
;



#endif /* XMIPPCORE_CORE_XMIPP_ARRAY_DIM_H_ */
