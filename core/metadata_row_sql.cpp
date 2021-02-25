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
    //Just initialize all pointers with NULL value
    memset(_objects, 0, MDL_LAST_LABEL * sizeof(size_t));
}

MDRowSql::MDRowSql(const MDRowSql &row) {
    _size = 0;
    //Just initialize all pointers with NULL value
    memset(_objects, 0, MDL_LAST_LABEL * sizeof(size_t));
    copy(row);
}

MDRowSql& MDRowSql::operator = (const MDRowSql &row) {
    copy(row);
    return *this;
}

MDRowSql::~MDRowSql() {
    MDObject **ptrObjectsLabel = &(_objects[0]);
    FOR_ALL_LABELS() {
        delete *ptrObjectsLabel;
        ++ptrObjectsLabel;
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
    FOR_ALL_LABELS() {
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
            res.push_back(label);
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

MDObject *MDRowSql::getObject(MDLabel label) const {
    return _objects[label];
}

bool MDRowSql::getValue(MDObject &object) const {
    int _label = object.label;
    if (_objects[_label] == nullptr)
        return false;
    object.copy(*(_objects[_label]));
    return true;
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
    FOR_ALL_LABELS() {
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
    memcpy(_order, row._order, sizeof(int)*_size);
}

MDObject* MDRowSql::iteratorValue(size_t i) const {
    return _objects[_order[i]];
}

std::ostream& operator << (std::ostream &out, const MDRowSql &row) {
    for (size_t i = 0; i < row._size; ++i) {
        row._objects[row._order[i]]->toStream(out);
        out << " ";
    }
    return out;
}

///////////////////////////////////////////////////////////////////////////////
// MDRowSqlConst

MDRowSqlConst::MDRowSqlConst() : _order() {
    std::fill(_order.begin(), _order.end(), -1);
}

MDRowSqlConst::MDRowSqlConst(const std::vector<MDObject> &values) : _values(values), _order() {
    std::fill(_order.begin(), _order.end(), -1);
    for (size_t i = 0; i < values.size(); i++)
        _order[values[i].label] = i;
}

MDRowSqlConst::MDRowSqlConst(const MDRowSqlConst &row)
    : _values(row._values), _order(row._order) {}

MDRowSqlConst& MDRowSqlConst::operator = (const MDRowSqlConst &row) {
    _values = row._values;
    _order = row._order;
    return *this;
}

bool MDRowSqlConst::empty() const {
    return _values.size() == 0;
}

int MDRowSqlConst::size() const {
    return _values.size();
}

bool MDRowSqlConst::containsLabel(MDLabel label) const {
    return _order[label] != -1;
}

std::vector<MDLabel> MDRowSqlConst::labels() const {
    std::vector<MDLabel> res;
    res.reserve(_values.size());
    for (const MDObject &obj : _values)
        res.push_back(obj.label);
    return res;
}

const MDObject *MDRowSqlConst::getObject(MDLabel label) const {
    return &_values[_order.at(label)];
}

bool MDRowSqlConst::getValue(MDObject &object) const {
    if (_order[object.label] == -1)
        return false;
    object.copy(_values[_order.at(object.label)]);
    return true;
}

const MDObject* MDRowSqlConst::iteratorValue(size_t i) const {
    return &_values[i];
}

std::ostream& operator << (std::ostream &out, const MDRowSqlConst &row) {
    for (size_t i = 0; i < row.size(); ++i) {
        row._values[i].toStream(out);
        out << " ";
    }
    return out;
}
