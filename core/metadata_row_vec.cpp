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
#include <cassert>
#include "metadata_row_vec.h"

MDRowVec::MDRowVec()
    : _col_to_label(nullptr), _no_columns(nullptr), _in_metadata(false) {
    _row = new std::vector<MDObject>();
    _label_to_col = new std::array<int, MDL_LAST_LABEL>();
    std::fill(_label_to_col->begin(), _label_to_col->end(), -1);
}

MDRowVec::MDRowVec(std::vector<MDObject>& row, size_t rowi, std::array<int, MDL_LAST_LABEL>& label_to_col,
                   std::vector<MDLabel>& col_to_label, size_t& no_columns)
    : _row(&row), _rowi(rowi), _label_to_col(&label_to_col), _col_to_label(&col_to_label),
      _no_columns(&no_columns), _in_metadata(true)
    {}

MDRowVec::MDRowVec(const std::vector<MDObject>& row, size_t rowi, const std::array<int, MDL_LAST_LABEL>& label_to_col,
                   const std::vector<MDLabel>& col_to_label, const size_t& no_columns)
    : _row(const_cast<std::vector<MDObject>*>(&row)),
       _rowi(rowi),
       _label_to_col(const_cast<std::array<int, MDL_LAST_LABEL>*>(&label_to_col)),
       _col_to_label(const_cast<std::vector<MDLabel>*>(&col_to_label)),
       _no_columns(const_cast<size_t*>(&no_columns)),
       _in_metadata(true)
       // This is very nasty hack. I was unable to solve situation nicely.
       // Creator of the object must pay close attention to create 'const MDRowVec' when instantiating from
       // const MetaData. This should be ok as crator is only in MetaDataVec.
    {}

MDRowVec::MDRowVec(const MDRowVec &other)
    : _rowi(other._rowi), _col_to_label(other._col_to_label), _no_columns(other._no_columns),
      _in_metadata(other._in_metadata)
    {
    if (_in_metadata) {
        _row = other._row;
        _label_to_col = other._label_to_col;
    } else {
        _row = new std::vector<MDObject>(*(other._row));
        _label_to_col = new std::array<int, MDL_LAST_LABEL>(*(other._label_to_col));
    }
}

MDRowVec MDRowVec::deepCopy(const MDRowVec& row) {
    MDRowVec newRow(row);
    if (newRow._in_metadata)
        newRow.detach();
    return newRow;
}

void MDRowVec::detach() {
    if (!this->_in_metadata)
        return;

    this->_in_metadata = false;
    this->_row = new std::vector<MDObject>(*(this->_row));
    this->_label_to_col = new std::array<int, MDL_LAST_LABEL>(*(this->_label_to_col));
    this->_col_to_label = nullptr;
    this->_no_columns = nullptr;
}

MDRowVec &MDRowVec::operator = (const MDRowVec &other) {
    if (this == &other)
        return *this;

    if (!_in_metadata) {
        delete _row;
        delete _label_to_col;
    }

    _rowi = other._rowi;
    _no_columns = other._no_columns;
    _in_metadata = other._in_metadata;
    _col_to_label = other._col_to_label;

    if (_in_metadata) {
        _row = other._row;
        _label_to_col = other._label_to_col;
    } else {
        _row = new std::vector<MDObject>(*(other._row));
        _label_to_col = new std::array<int, MDL_LAST_LABEL>(*(other._label_to_col));
    }

    return *this;
}

MDRow& MDRowVec::operator = (const MDRow& row) {
    *this = dynamic_cast<const MDRowVec&>(row);
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
    if (!_in_metadata)
        std::fill(_label_to_col->begin(), _label_to_col->end(), -1);
}

bool MDRowVec::containsLabel(MDLabel label) const {
    return (*_label_to_col)[label] >= 0;
}

std::vector<MDLabel> MDRowVec::labels() const {
    std::vector<MDLabel> res;
    res.reserve(_row->size());
    for (const auto mdObj: *_row)
        res.emplace_back(mdObj.label);
    return res;
}

void MDRowVec::addLabel(MDLabel label) {
    // Warning: not adding to all rows!
    if ((*_label_to_col)[label] < 0)
        newCol(label);
}

size_t MDRowVec::newCol(const MDLabel label) {
    size_t i;
    if (_no_columns != nullptr)
        i = (*_no_columns)++;
    else
        i = _row->size();

    if (_col_to_label != nullptr) {
        while (_row->size() < i) {
            size_t j = _row->size();
            _row->emplace_back(MDObject((*_col_to_label).at(j)));
        }
        if (i < (*_col_to_label).size()) {
            (*_col_to_label).at(i) = label;
        } else {
            assert(i == (*_col_to_label).size());
            (*_col_to_label).push_back(label);
        }
    }

    assert(_row->size() == i);
    (*_label_to_col).at(label) = i;
    _row->emplace_back(MDObject(label));
    return i;
}

MDObject *MDRowVec::getObject(MDLabel label) {
    if ((*_label_to_col)[label] < 0)
        return nullptr;
    if ((*_label_to_col)[label] >= static_cast<int>(_row->size()))
        return nullptr;
    return &_row->at((*_label_to_col)[label]);
}

const MDObject *MDRowVec::getObject(MDLabel label) const {
    if ((*_label_to_col)[label] < 0)
        return nullptr;
    if ((*_label_to_col)[label] >= static_cast<int>(_row->size()))
        return nullptr;
    return &_row->at((*_label_to_col)[label]);
}

void MDRowVec::setValue(const MDObject &object) {
    MDLabel _label = object.label;
    int coli = (*_label_to_col)[_label];
    if (coli < 0) {
        size_t i = newCol(_label);
        (*_row)[i] = object;
    } else {
        while (_row->size() <= coli) {
            size_t j = _row->size();
            _row->emplace_back(MDObject((*_col_to_label)[j]));
        }

        (*_row)[(*_label_to_col)[_label]] = object;
    }
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
