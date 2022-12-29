/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#ifndef CORE_RWHDF5_H_
#define CORE_RWHDF5_H_


///@defgroup HDF5 HDF5 File format
///@ingroup ImageFormats
//@{


/** Determine datatype of a HDF5 dataset.
  * @ingroup TIFF
  */
DataType datatypeH5(hid_t dataset);

hid_t H5Datatype(DataType datatype);


/** Read Images from HDF5 container files.
  */
int readHDF5(size_t select_img);

/** Write Images to HDF5 container files.
  */
int writeHDF5(size_t select_img, bool isStack=false, int mode=WRITE_OVERWRITE, String bitDepth="", CastWriteMode castMode = CW_CAST);

//@}
#endif /* RWHDF5_H_ */
