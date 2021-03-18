/***************************************************************************
 *
 * Authors:    Jan Horacek (xhorace4@fi.muni.cz)
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

#include <algorithm>
#include <sstream>
#include "metadata_row_vec.h"

MDRowVec::MDRowVec()
    : _in_metadata(false) {
    _row = new std::vector<MDObject>();
    _label_to_col = new std::array<int, MDL_LAST_LABEL>();
}

MDRowVec::MDRowVec(std::vector<MDObject>& row, size_t rowi, std::array<int, MDL_LAST_LABEL>& label_to_col)
    : _row(&row), _rowi(rowi), _label_to_col(&label_to_col), _in_metadata(true)
    {}

MDRowVec::MDRowVec(const std::vector<MDObject>& row, size_t rowi, const std::array<int, MDL_LAST_LABEL>& label_to_col)
    : _row(const_cast<std::vector<MDObject>*>(&row)),
       _rowi(rowi),
       _label_to_col(const_cast<std::array<int, MDL_LAST_LABEL>*>(&label_to_col)), _in_metadata(true)
       // This is very nasty hack. I was unable to solve situation nicely.
       // Creator of the object must pay close attention to create 'const MDRowVec' when instantiating from
       // const MetaData. This should be ok as crator is only in MetaDataVec.
    {}

MDRowVec::MDRowVec(const MDRowVec &other)
    : _rowi(other._rowi), _in_metadata(other._in_metadata)
    {
    if (_in_metadata) {
        _row = other._row;
        _label_to_col = other._label_to_col;
    } else {
        _row = new std::vector<MDObject>(*(other._row));
        _label_to_col = new std::array<int, MDL_LAST_LABEL>(*(other._label_to_col));
    }
}

MDRowVec &MDRowVec::operator = (const MDRowVec &other) {
    _rowi = other._rowi;
    _in_metadata = other._in_metadata;

    if (_in_metadata) {
        _row = other._row;
        _label_to_col = other._label_to_col;
    } else {
        _row = new std::vector<MDObject>(*(other._row));
        _label_to_col = new std::array<int, MDL_LAST_LABEL>(*(other._label_to_col));
    }

    return *this;
}

MDRowVec::~MDRowVec() {
    if (!_in_metadata) {
        delete _row;
        delete _label_to_col;
    }
}

bool MDRowVec::empty() const {
    return _row->size() == 0;
}

int MDRowVec::size() const {
    return _row->size();
}

void MDRowVec::clear() {
    _row->clear();
}

bool MDRowVec::containsLabel(MDLabel label) const {
    return (*_label_to_col)[label] >= 0;
}

std::vector<MDLabel> MDRowVec::labels() const {
    std::vector<MDLabel> res;
    res.reserve(_row->size());
    for (const auto mdObj: *_row)
        res.push_back(mdObj.label);
    return res;
}

void MDRowVec::addLabel(MDLabel label) {
    // Warning: not adding to all rows!
    if ((*_label_to_col)[label] < 0) {
        (*_label_to_col)[label] = _row->size();
        _row->push_back(MDObject(label));
    }
}

MDObject *MDRowVec::getObject(MDLabel label) {
    return &_row->at((*_label_to_col)[label]);
}

const MDObject *MDRowVec::getObject(MDLabel label) const {
    return &_row->at((*_label_to_col)[label]);
}

bool MDRowVec::getValue(MDObject &object) const {
    MDLabel _label = object.label;
    if ((*_label_to_col)[_label] < 0)
        return false;
    object.copy(_row->at((*_label_to_col)[_label]));
    return true;
}

void MDRowVec::setValue(const MDObject &object) {
    MDLabel _label = object.label;
    if ((*_label_to_col)[_label] < 0) {
        (*_label_to_col)[_label] = _row->size();
        _row->push_back(object);
    } else
        (*_row)[(*_label_to_col)[_label]] = object;
}

MDObject* MDRowVec::iteratorValue(size_t i) {
    return &(*_row)[i];
}

const MDObject* MDRowVec::iteratorValue(size_t i) const {
    return &(*_row)[i];
}

std::ostream& operator << (std::ostream &out, const MDRowVec &row) {
    for (int i = 0; i < row.size(); ++i) {
        (*(row._row))[i].toStream(out);
        out << " ";
    }
    return out;
}

bool MDRowVec::inMetadata() const { return _in_metadata; }
