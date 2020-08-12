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

#ifndef XMIPPCORE_CORE_XMIPP_WRITE_MODE_H_
#define XMIPPCORE_CORE_XMIPP_WRITE_MODE_H_

/** Write mode
 * This class defines the writing behavior.
 */
typedef enum
{
    WRITE_READONLY,   //only can read the file
    WRITE_OVERWRITE, //forget about the old file and overwrite it
    WRITE_REPLACE,   //replace a particular object by another
    WRITE_APPEND,    //append and object at the end of a stack, so far can not append stacks
    WRITE_LAST_LABEL                       // **** NOTE ****: Do keep this label always at the end
    // it is here for looping purposes
} WriteMode;



#endif /* XMIPPCORE_CORE_XMIPP_WRITE_MODE_H_ */
