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

#ifndef CORE_XMIPP_IMAGE_FHANDLER_H_
#define CORE_XMIPP_IMAGE_FHANDLER_H_

#include <hdf5.h>
#include <tiffio.h>

/** Open File struct
 * This struct is used to share the File handlers with Image Collection class
 */
struct ImageFHandler
{
    FILE*     fimg;       // Image File handler
    FILE*     fhed;       // Image File header handler
    TIFF*     tif;        // TIFF Image file handler
    hid_t     fhdf5;   // HDF5 File handler
    FileName  fileName;   // Image file name
    FileName  headName;   // Header file name
    FileName  ext_name;   // Filename extension
    bool     exist;       // Shows if the file exists. Equal 0 means file does not exist or not stack.
    int        mode;   // Opening mode behavior
};


#endif /* CORE_XMIPP_IMAGE_FHANDLER_H_ */
