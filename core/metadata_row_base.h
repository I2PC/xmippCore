/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Jan Horacek (xhorace4@fi.muni.cz)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 * Institute of Computer Science MUNI
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

#ifndef CORE_METADATA_ROW_BASE_H
#define CORE_METADATA_ROW_BASE_H

#include "metadata_label.h"
#include "metadata_object.h"
#include "choose.h"
#include <stdexcept>

/** Get value */
//    template <typename T >
//    void getValueOrAbort(MDLabel label, T &d) const
//    {
//        if (!getValue(label, d))
//         REPORT_ERROR(ERR_ARG_MISSING,(String)"Cannot find label: " + MDL::label2Str(label) );
//        //formatString("%d",label) );
//    }
//weite function as macro since MDL::label2Str is not availale at class compilation time
#define rowGetValueOrAbort(__row,__label, __d)\
    if (!__row.getValue(__label, __d))\
     REPORT_ERROR(ERR_ARG_MISSING,(String)"Cannot find label: " + MDL::label2Str(__label) );


/** Common API of all metadata rows (abstract class).
 * Classes like MDRowVec & MDRowSql implement this API.
 */
class MDRow {
public:
    /* Row could be attached to metadata (contains pointers to MetaData) or detached
     * from metadata (contains data itself). To detach row, call this method.
     */
    virtual void detach() {}

    virtual bool empty() const = 0;
    virtual int size() const = 0; /** Return number of labels present */
    virtual void clear() = 0;

    virtual MDRow& operator = (const MDRow&) = 0;
    virtual ~MDRow() = default;

    virtual size_t id() const {
        return getObject(MDL_OBJID)->getValue2(size_t());
    }

    virtual bool containsLabel(MDLabel label) const = 0;
    virtual std::vector<MDLabel> labels() const = 0;
    virtual void addLabel(MDLabel label) = 0;

    /** Reset the values of the labels related to
     *  geometry to their default values
     */
    virtual void resetGeo(bool addLabels = true) {
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

    virtual MDObject* getObject(MDLabel label) = 0;
    virtual const MDObject* getObject(MDLabel label) const = 0;

    template <typename T>
    T& getValue(MDLabel label) {
        MDObject* ptr = getObject(label);
        if (ptr == nullptr)
            throw std::logic_error("Object does not exist!");
        return ptr->getValue2(T());
    }

    template <typename T>
    const T& getValue(MDLabel label) const {
        const MDObject* ptr = getObject(label);
        if (ptr == nullptr)
            throw std::logic_error("Object does not exist!");
        return ptr->getValue2(T());
    }

    template <typename T>
    bool getValue(MDLabel label, T &d) const { // FIXME: deprecated
        auto *obj = getObject(label);
        if (obj == nullptr)
            return false;
        d = obj->getValue2(T());
        return true;
    }

    bool getValue(MDObject &object) const { // FIXME: deprecated
        const MDObject* ptr = this->getObject(object.label);
        if (ptr != nullptr)
            object = *ptr;
        return (ptr != nullptr);
    }

    template <typename T>
    const T& getValueOrDefault(MDLabel label, const T& def) const {
        const MDObject* ptr = getObject(label);
        if (ptr == nullptr)
            return def;
        return ptr->getValue2(T());
    }

    template <typename T>
    T& getValueOrDefault(MDLabel label, const T& def) {
        const MDObject* ptr = getObject(label);
        if (ptr == nullptr)
            return def;
        return ptr->getValue2(T());
    }

    template <typename T, typename T1>
    void getValueOrDefault(MDLabel label, T &d, T1 def) const { // FIXME: deprecated
        if (!getValue(label, d))
            d = (T) def;
    }

    template <typename T>
    void setValue(MDLabel label, const T &d, bool addLabel = true) {
        if (this->containsLabel(label) || addLabel)
            this->setValue(MDObject(label, d));
    }

    virtual void setValue(const MDObject &object) = 0;

    virtual void setValueFromStr(MDLabel label, const String &value) {
        MDObject mdValue(label);
        mdValue.fromString(value);
        setValue(mdValue);
    }

    friend std::ostream& operator << (std::ostream &out, const MDRow &row);

    template <bool IsConst>
    class iterator_ptr {
        size_t i;
        const MDRow* row;

    public:
        iterator_ptr(size_t i, const MDRow& row) : i(i), row(&row) {}
        iterator_ptr(iterator_ptr const& right) : i(right.i), row(right.row) {}
        iterator_ptr& operator=(iterator_ptr const& right) {
            i = right.i;
            row = right.row;
            return *this;
        }
        iterator_ptr& operator++() {
            ++i;
            return *this;
        }
        typename TypeHelpers::choose<IsConst, const MDObject*, MDObject*>::type operator*() const {
            return row->iteratorValue(i);
        }
        bool operator==(const iterator_ptr& other) const { return other.i == this->i; }
        bool operator!=(const iterator_ptr& other) const { return !(*this == other); }
    };

    using iterator = iterator_ptr<false>;
    using const_iterator = iterator_ptr<true>;

    iterator begin() { return iterator_ptr<false>(0, *this); }
    iterator end() { return iterator_ptr<false>(this->size(), *this); }

    const_iterator begin() const { return iterator_ptr<true>(0, *this); }
    const_iterator end() const { return iterator_ptr<true>(this->size(), *this); }

private:

    virtual MDObject* iteratorValue(size_t i) = 0;
    virtual const MDObject* iteratorValue(size_t i) const = 0;
};

#endif
