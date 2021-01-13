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

/*bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label) {
    std::vector<MDLabel>::const_iterator location;
    location = std::find(labelsVector.begin(), labelsVector.end(), label);
    return (location != labelsVector.end());
}*/

MDRowVec::MDRowVec(const MDRowVec &other)
    : _row(other._row), _rowi(other._rowi), _label_to_col(other._label_to_col)
    {}

MDRowVec &MDRowVec::operator = (const MDRowVec &other) {
    _row = other._row;
    _rowi = other._rowi;
    _label_to_col = other._label_to_col;
}

bool MDRowVec::empty() const {
    return _row.size() == 0;
}

int MDRowVec::size() const {
    return _row.size();
}

void MDRowVec::clear() {
    _row.clear();
}

bool MDRowVec::containsLabel(MDLabel label) const {
    return _label_to_col[label] >= 0;
}

std::vector<MDLabel> MDRowVec::labels() const {
    std::vector<MDLabel> res;
    res.reserve(_row.size());
    for (const auto mdObj: _row)
        res.push_back(mdObj.label);
    return res;
}

void MDRowVec::addLabel(MDLabel label) {
    // Warning: not adding to all rows!
    if (_label_to_col[label] < 0) {
        _label_to_col[label] = _row.size();
        _row.push_back(MDObject(label));
    }
}

MDObject *MDRowVec::getObject(MDLabel label) const {
    return &_row.at(_label_to_col[label]);
}

bool MDRowVec::getValue(MDObject &object) const {
    MDLabel _label = object.label;
    if (_label_to_col[_label] < 0)
        return false;
    object.copy(_row.at(_label_to_col[_label]));
    return true;
}

void MDRowVec::setValue(const MDObject &object) {
    MDLabel _label = object.label;
    if (_label_to_col[_label] < 0) {
        _label_to_col[_label] = _row.size();
        _row.push_back(object);
    } else
        _row[_label_to_col[_label]] = object;
}

MDObject* MDRowVec::iteratorValue(size_t i) const {
    return &_row[i];
}

std::ostream& operator << (std::ostream &out, const MDRowVec &row) {
    for (int i = 0; i < row.size(); ++i) {
        row._row[i].toStream(out);
        out << " ";
    }
    return out;
}
