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

#ifndef XMIPPCORE_CORE_TRANSFORMATIONS_DEFINES_H_
#define XMIPPCORE_CORE_TRANSFORMATIONS_DEFINES_H_


namespace xmipp_transformation
{
	enum XmippInterpolation {NEAREST=0, LINEAR=1, BSPLINE2=2, BSPLINE3=3, BSPLINE4=4}; // Interpolation type
	static const bool IS_INV=true;
	static const bool IS_NOT_INV=false;
	static const bool DONT_WRAP=false;
	static const bool WRAP=true;
}

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-50
#endif

#endif /* XMIPPCORE_CORE_TRANSFORMATIONS_DEFINES_H_ */
