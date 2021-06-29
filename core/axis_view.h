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

#ifndef CORE_AXIS_VIEW_H_
#define CORE_AXIS_VIEW_H_

/**
 * Possible views for 3D MuldimArray
 */
typedef enum
{
    VIEW_Z_NEG,    // Front view (Z negative)
    VIEW_Z_POS,    //  Z positve
    VIEW_Y_NEG,    // Align -Y axis to Z axis, rotating 90 degrees around X axis");
    VIEW_Y_POS, // Align Y axis to Z axis, rotating -90 degrees around X axis");
    VIEW_X_NEG,   // Align -X axis to Z axis, rotating -90 degrees around Y axis");
    VIEW_X_POS   // Align X axis to Z axis, rotating 90 degrees around Y axis");
} AxisView;


#endif /* CORE_AXIS_VIEW_H_ */
