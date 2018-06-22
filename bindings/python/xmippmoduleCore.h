/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Roberto Marabini       (roberto@cnb.csic.es)
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

#ifndef _XMIPPMODULECORE_H
#define _XMIPPMODULECORE_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION // this code is NumPy 1.8 compliant (i.e. we don't need API deprecated in 1.7)

#include <Python.h>

#include "../../core/metadata_extension.h"
#include "../../core/xmipp_image_generic.h"
#include "../../core/xmipp_image_extension.h"
#include "../../core/xmipp_color.h"
#include "../../core/symmetries.h"

#include "python_filename.h"
#include "python_image.h"
#include "python_program.h"
#include "python_metadata.h"
#include "python_symmetry.h"

extern PyObject * PyXmippError;

#define SymList_Check(v) (((v)->ob_type == &SymListType))
#define SymList_Value(v)  ((*((SymListObject*)(v))->symlist))


/***************************************************************/
/*                            Global methods                   */
/***************************************************************/
PyObject *
xmipp_str2Label(PyObject *obj, PyObject *args);

PyObject *
xmipp_label2Str(PyObject *obj, PyObject *args);

PyObject *
xmipp_colorStr(PyObject *obj, PyObject *args);

PyObject *
xmipp_labelType(PyObject *obj, PyObject *args);

PyObject *
xmipp_labelHasTag(PyObject *obj, PyObject *args);

PyObject *
xmipp_labelIsImage(PyObject *obj, PyObject *args);

/* isInStack */
PyObject *
xmipp_isValidLabel(PyObject *obj, PyObject *args, PyObject *kwargs);

/* createEmptyFile */
PyObject *
xmipp_createEmptyFile(PyObject *obj, PyObject *args, PyObject *kwargs);

/* getImageSize */
PyObject *
xmipp_getImageSize(PyObject *obj, PyObject *args, PyObject *kwargs);

/* ImgSize (from metadata filename)*/
PyObject *
xmipp_ImgSize(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_CheckImageFileSize(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_CheckImageCorners(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_ImgCompare(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_compareTwoFiles(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_compareTwoImageTolerance(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_bsoftRemoveLoopBlock(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_bsoftRestoreLoopBlock(PyObject *obj, PyObject *args, PyObject *kwargs);


/***************************************************************/
/*                   Some specific utility functions           */
/***************************************************************/
/* readMetaDataWithTwoPossibleImages */
PyObject *
xmipp_readMetaDataWithTwoPossibleImages(PyObject *obj, PyObject *args,
                                        PyObject *kwargs);

/* substituteOriginalImages */
PyObject *
xmipp_substituteOriginalImages(PyObject *obj, PyObject *args, PyObject *kwargs);

bool validateInputImageString(PyObject * pyImage, PyObject *pyStrFn, FileName &fn);

PyObject *
xmipp_compareTwoMetadataFiles(PyObject *obj, PyObject *args, PyObject *kwargs);

/* dump metadatas to database*/
PyObject *
xmipp_dumpToFile(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_Euler_angles2matrix(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
xmipp_Euler_direction(PyObject *obj, PyObject *args, PyObject *kwargs);

PyObject *
MetaData_activateMathExtensions(PyObject *obj, PyObject *args, PyObject *kwargs);

void addIntConstant(PyObject * dict, const char * name, const long &value);
void addLabels(PyObject * dict);

#endif
