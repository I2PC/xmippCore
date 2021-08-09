/***************************************************************************
 *
 * Authors:     Jan Horacek (xhorace4@fi.muni.cz)
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

/** Definitions of API of MetaData iterators (abstract classes). */

#ifndef IT_BASE_METADATA_H
#define IT_BASE_METADATA_H

#include <memory>
#include "metadata_row_base.h"
#include "choose.h"

/** Iterates over metadata rows */
template <bool IsConst>
struct MDBaseRowIterator {
    virtual ~MDBaseRowIterator() {}
    virtual std::unique_ptr<MDBaseRowIterator> clone() = 0;
    virtual void increment() = 0;
    virtual bool operator==(const MDBaseRowIterator& other) const = 0;
    virtual bool operator!=(const MDBaseRowIterator& other) const { return !(*this == other); }
    virtual typename TypeHelpers::choose<IsConst, const MDRow&, MDRow&>::type operator*() = 0;
};

/** Iterates over metadata ids */
template <bool IsConst>
struct MDBaseIdIterator {
    virtual ~MDBaseIdIterator() {}
    virtual std::unique_ptr<MDBaseIdIterator> clone() = 0;
    virtual void increment() = 0;
    virtual bool operator==(const MDBaseIdIterator& other) const = 0;
    virtual bool operator!=(const MDBaseIdIterator& other) const { return !(*this == other); }
    virtual size_t operator*() = 0;
};

#endif
