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

/** Class for holding an entire row of posible MDObject */
class MDRow
{
public:
    //Reserve space for the maximum different labels
    //this will allows constant access to each object indexing by labels
    MDObject * objects[MDL_LAST_LABEL];
    MDLabel order[MDL_LAST_LABEL];
    int _size; //Number of active labels

public:
    /** Empty constructor */
    MDRow();
    /** Copy constructor */
    MDRow(const MDRow & row);

    /** Assignment */
    MDRow& operator = (const MDRow &row);

    /** Destructor */
    ~MDRow();
    /** True if this row contains this label */
    inline bool containsLabel(MDLabel label) const {
        return objects[label] != NULL;
    }

    /** Returns all labels defined for this row */
    std::vector<MDLabel> getLabels() const;

    /** Add a new label */
    void addLabel(MDLabel label);

    /** Clear elements of the row */
    void clear();
    /** Return number of labels present */
    int size() const;
    /** Function to test whether is empty */
    bool empty() const;

    /** Reset the values of the labels related to
     *  geometry to their default values
     */
    void resetGeo(bool addLabels = true);

    /** Get object */
    MDObject * getObject(MDLabel label) const;

    /** Get value */
    template <typename T>
    bool getValue(MDLabel label, T &d) const
    {
        if (objects[label] == NULL)
            return false;

        objects[label]->getValue(d);
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

    /** Get value */
    template <typename T, typename T1>
    void getValueOrDefault(MDLabel label, T &d, T1 def) const
    {
        if (!getValue(label, d))
            d = (T) def;
    }

    bool getValue(MDObject &object) const;

    /** Set value
     *
     * @param label    Metadata label to be set
     * @param d        Value of the label
     * @param addLabel Add label if is not already contained
     */
    template <typename T>
    void setValue(MDLabel label, const T &d, bool addLabel = true)
    {
        if (objects[label] != NULL)
            objects[label]->setValue(d);
        else if (addLabel)
        {
            objects[label] = new MDObject(label, d);
            order[_size] = label;
            ++_size;
        }
    }

    void setValue(const MDObject &object);

    void setValueFromStr(MDLabel label, const String &value);

    /** Show */
    friend std::ostream& operator << (std::ostream &out, const MDRow &row);

private:
    void copy(const MDRow &row);
};

#endif
