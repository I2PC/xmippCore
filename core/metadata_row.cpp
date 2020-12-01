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
#include "metadata_row_base.h"

bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label)
{
    std::vector<MDLabel>::const_iterator location;
    location = std::find(labelsVector.begin(), labelsVector.end(), label);

    return (location != labelsVector.end());
}

void MDRow::clear()
{
    _size = 0;
    //Just initialize all pointers with NULL value
    FOR_ALL_LABELS()
    {
        delete objects[_label];
        objects[_label] = NULL;
    }
}

bool MDRow::empty() const
{
    return _size == 0;
}

void MDRow::resetGeo(bool addLabels)
{
    setValue(MDL_ORIGIN_X,  0., addLabels);
    setValue(MDL_ORIGIN_Y,  0., addLabels);
    setValue(MDL_ORIGIN_Z,  0., addLabels);
    setValue(MDL_SHIFT_X,   0., addLabels);
    setValue(MDL_SHIFT_Y,   0., addLabels);
    setValue(MDL_SHIFT_Z,   0., addLabels);
    setValue(MDL_ANGLE_ROT, 0., addLabels);
    setValue(MDL_ANGLE_TILT,0., addLabels);
    setValue(MDL_ANGLE_PSI, 0., addLabels);
    setValue(MDL_WEIGHT,   1., addLabels);
    setValue(MDL_FLIP,     false, addLabels);
    setValue(MDL_SCALE,    1., addLabels);
}

int MDRow::size() const
{
    return _size;
}

std::vector<MDLabel> MDRow::getLabels() const {
    std::vector<MDLabel> res;
    res.reserve(_size);
    for (int i = 0; i < _size; ++i){
        const MDLabel &label = order[i];
        if (containsLabel(label)) {
            res.push_back(label);
        }
    }
    return res;
}

void MDRow::addLabel(MDLabel label)
{
    if (objects[label] == NULL)
    {
        objects[label] = new MDObject(label);
        order[_size] = label;
        ++_size;
    }
}

MDObject * MDRow::getObject(MDLabel label) const
{
    return objects[label];
}

/** Get value */
bool MDRow::getValue(MDObject &object) const
{
    int _label = object.label;
    if (objects[_label] == NULL)
        return false;
    object.copy(*(objects[_label]));
    return true;
}

/** Useful macro for copy values */
/** Set value */
void MDRow::setValue(const MDObject &object)
{
    int _label = object.label;
    if (objects[_label] == NULL)
    {
        objects[_label] = new MDObject(object);
        order[_size] = object.label;
        ++_size;
    }
    else
        objects[_label]->copy(object);
}

void MDRow::setValueFromStr(MDLabel label, const String &value)
{
    MDObject mdValue(label);
    mdValue.fromString(value);
    setValue(mdValue);
}

MDRow::~MDRow()
{
    MDObject ** ptrObjectsLabel=&(objects[0]);
    FOR_ALL_LABELS()
    {
        delete *ptrObjectsLabel;
        ++ptrObjectsLabel;
    }
}

MDRow::MDRow(const MDRow & row)
{
    _size = 0;
    //Just initialize all pointers with NULL value
    memset(objects, 0, MDL_LAST_LABEL * sizeof(size_t));
    copy(row);
}

MDRow::MDRow()
{
    _size = 0;
    //Just initialize all pointers with NULL value
    memset(objects, 0, MDL_LAST_LABEL * sizeof(size_t));
}

MDRow& MDRow::operator = (const MDRow &row)
{
    copy(row);
    return *this;
}

void MDRow::copy(const MDRow &row)
{
    //Copy existing MDObjects from row
    //and delete unexisting ones
    _size = row._size;
    MDObject ** ptrObjectsLabel=&(objects[0]);
    MDObject * const * ptrRowObjectsLabel=&(row.objects[0]);
    FOR_ALL_LABELS()
    {
        if (*ptrRowObjectsLabel == NULL)
        {
            delete *ptrObjectsLabel;
            *ptrObjectsLabel = NULL;
        }
        else
        {
            if (*ptrObjectsLabel == NULL)
                *ptrObjectsLabel = new MDObject(*(*ptrRowObjectsLabel));
            else
                (*ptrObjectsLabel)->copy(*(*ptrRowObjectsLabel));
        }
        ++ptrObjectsLabel;
        ++ptrRowObjectsLabel;
    }
    //copy the order of labels
    memcpy(order, row.order, sizeof(int)*_size);
}

std::ostream& operator << (std::ostream &out, const MDRow &row)
{
    for (int i = 0; i < row._size; ++i)
    {
        row.objects[row.order[i]]->toStream(out);
        out << " ";
    }
    return out;
}
