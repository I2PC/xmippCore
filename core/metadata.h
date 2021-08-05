/***************************************************************************
 *
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Jan Horacek (xhorace4@fi.muni.cz)
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

/** MetaData class implements a way to store metadata in xmipp programs.
 * - MetaData is represented as as a database-like table.
 * - Each intersection of row and column contains MDObject instance.
 * - Each row has its label MDLabel.
 * - Columns of all rows are same.
 * - Metadata could be loaded from file, stored to file, iterated over, rows
 *   can be added, removed, changed etc.
 *
 * There was an original database (MetaDataDb) implementation of MetaData in
 * this file till 2021, when MetaData were split into MetaDataDb & MetaDataVec
 * with aim to achive higner speeds via saving metadata in std::vectors instead
 * of sql database.
 *
 * This file is present for backward-compatibility only and as general
 * description of metadata.

 * Current MetaData implementation:
 *  1. metadata_base.(h|cpp): common MetaData API definition
 *  2. metadata_vec.(h|cpp): vector MetaData implementation
 *  3. metadata_db.(h|cpp): old datababse MetaData implementation
 */

#include "metadata_db.h"
