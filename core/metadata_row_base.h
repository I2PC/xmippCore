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

struct MDRowIterator {
};

/** Abstract class (API) for holding an entire row of posible MDObject */
class MDRow {
public:
    using iterator = MDRowIterator;

    virtual bool empty() const = 0;
    virtual int size() const = 0; /** Return number of labels present */
    virtual void clear() = 0;

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

    virtual MDObject *getObject(MDLabel label) const = 0;

    template <typename T>
    bool getValue(MDLabel label, T &d) const {
        if (getObject(label) == nullptr)
            return false;
        getObject(label)->getValue(d);
        return true;
    }

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

    template <typename T, typename T1>
    void getValueOrDefault(MDLabel label, T &d, T1 def) const {
        if (!getValue(label, d))
            d = (T) def;
    }

    virtual bool getValue(MDObject &object) const = 0;

    /** Set value
     *
     * @param label    Metadata label to be set
     * @param d        Value of the label
     * @param addLabel Add label if is not already contained
     */
    template <typename T>
    void setValue(MDLabel label, const T &d, bool addLabel = true) {
        MDObject *obj = getObject(label);
        if (obj == nullptr && addLabel) {
            this->addLabel(label);
            obj = getObject(label);
        }

        if (obj != nullptr)
            obj->setValue(d);
    }

    virtual void setValue(const MDObject &object) = 0;

    virtual void setValueFromStr(MDLabel label, const String &value) {
        MDObject mdValue(label);
        mdValue.fromString(value);
        setValue(mdValue);
    }

    friend std::ostream& operator << (std::ostream &out, const MDRow &row);

private:

    virtual const MDLabel& order(size_t i) const = 0;
};

#endif
