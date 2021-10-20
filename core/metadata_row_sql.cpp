/***************************************************************************
 *
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *             Jan Horacek (xhorace4@fi.muni.cz)
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
#include "metadata_row_sql.h"

MDRowSql::MDRowSql() {
    _size = 0;
    _objects.fill(nullptr);
}

MDRowSql::MDRowSql(const MDRowSql &row) {
    _size = 0;
    _objects.fill(nullptr);
    copy(row);
}

MDRowSql& MDRowSql::operator = (const MDRowSql &row) {
    copy(row);
    return *this;
}

MDRow& MDRowSql::operator = (const MDRow& row) {
    *this = dynamic_cast<const MDRowSql&>(row);
    return *this;
}

MDRowSql::MDRowSql(const std::vector<MDObject> &values) {
    _objects.fill(nullptr);
    for (size_t i = 0; i < values.size(); i++) {
        _objects[values[i].label] = new MDObject(values[i]);
        _order[i] = values[i].label;
    }
}

MDRowSql::~MDRowSql() {
    for (int _label = MDL_FIRST_LABEL; _label < MDL_LAST_LABEL; ++_label) {
        if (_objects[_label] != nullptr)
            delete _objects[_label];
    }
}

bool MDRowSql::empty() const {
    return _size == 0;
}

int MDRowSql::size() const {
    return _size;
}

void MDRowSql::clear() {
    _size = 0;
    for (int _label = MDL_FIRST_LABEL; _label < MDL_LAST_LABEL; ++_label) {
        if (_objects[_label] != nullptr)
            delete _objects[_label];
        _objects[_label] = nullptr;
    }
}

bool MDRowSql::containsLabel(MDLabel label) const {
    return _objects[label] != nullptr;
}

std::vector<MDLabel> MDRowSql::labels() const {
    std::vector<MDLabel> res;
    res.reserve(_size);
    for (size_t i = 0; i < _size; ++i) {
        const MDLabel &label = _order[i];
        if (containsLabel(label))
            res.emplace_back(label);
    }
    return res;
}

void MDRowSql::addLabel(MDLabel label) {
    if (_objects[label] == nullptr) {
        _objects[label] = new MDObject(label);
        _order[_size] = label;
        ++_size;
    }
}

const MDObject *MDRowSql::getObject(MDLabel label) const {
    return _objects[label];
}

MDObject *MDRowSql::getObject(MDLabel label) {
    return _objects[label];
}

void MDRowSql::setValue(const MDObject &object) {
    int _label = object.label;
    if (_objects[_label] == nullptr) {
        _objects[_label] = new MDObject(object);
        _order[_size] = object.label;
        ++_size;
    } else
        _objects[_label]->copy(object);
}

void MDRowSql::copy(const MDRowSql &row) {
    //Copy existing MDObjects from row
    //and delete unexisting ones
    _size = row._size;
    MDObject ** ptrObjectsLabel=&(_objects[0]);
    MDObject * const * ptrRowObjectsLabel=&(row._objects[0]);
    for (int _label = MDL_FIRST_LABEL; _label < MDL_LAST_LABEL; ++_label) {
        if (*ptrRowObjectsLabel == nullptr) {
            delete *ptrObjectsLabel;
            *ptrObjectsLabel = nullptr;
        } else {
            if (*ptrObjectsLabel == nullptr)
                *ptrObjectsLabel = new MDObject(*(*ptrRowObjectsLabel));
            else
                (*ptrObjectsLabel)->copy(*(*ptrRowObjectsLabel));
        }
        ++ptrObjectsLabel;
        ++ptrRowObjectsLabel;
    }
    //copy the order of labels
    _order = row._order;
}

MDObject* MDRowSql::iteratorValue(size_t i) {
    return _objects[_order[i]];
}

const MDObject* MDRowSql::iteratorValue(size_t i) const {
    return _objects[_order[i]];
}

std::ostream& operator << (std::ostream &out, const MDRowSql &row) {
    for (size_t i = 0; i < row._size; ++i) {
        row._objects[row._order[i]]->toStream(out);
        out << " ";
    }
    return out;
}
