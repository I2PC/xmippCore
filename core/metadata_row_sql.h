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

#ifndef CORE_METADATA_ROW_SQL_H
#define CORE_METADATA_ROW_SQL_H

#include "metadata_row_base.h"

#include <array>

/** Class for holding an entire row of MDObject in MetaDataDb. */
class MDRowSql : public MDRow {
private:
    // Reserve space for the maximum different labels
    // this will allow constant access to each object indexing by labels
    std::array<MDObject*, MDL_LAST_LABEL> _objects; // label to object; nullptr if none
    std::array<MDLabel, MDL_LAST_LABEL> _order; // index to label (_order[0] = label for column 0)
    size_t _size; // Number of active labels
    size_t _id;

    void copy(const MDRowSql &row);
    MDObject* iteratorValue(size_t i) override;
    const MDObject* iteratorValue(size_t i) const override;

public:
    MDRowSql();
    ~MDRowSql();
    MDRowSql(const MDRowSql &row);
    MDRowSql(const std::vector<MDObject> &values);
    MDRow& operator = (const MDRow&) override;
    MDRowSql& operator = (const MDRowSql &row);

    size_t id() const override { return this->_id; }
    void set_id(size_t id) { this->_id = id; }

    bool empty() const override;
    int size() const override;
    void clear() override;

    bool containsLabel(MDLabel label) const override;
    std::vector<MDLabel> labels() const override;
    void addLabel(MDLabel label) override;

    MDObject *getObject(MDLabel label) override;
    const MDObject *getObject(MDLabel label) const override;

    void setValue(const MDObject &object) override;

    friend std::ostream& operator << (std::ostream &out, const MDRowSql &row);

    // Templated functions from base class must be retemplated

    template <typename T>
    T& getValue(MDLabel label) { return MDRow::getValue<T>(label); }

    template <typename T>
    const T& getValue(MDLabel label) const { return MDRow::getValue<T>(label); }

    template <typename T>
    bool getValue(MDLabel label, T &d) const { // FIXME: deprecated
        return MDRow::getValue<T>(label, d);
    }

    template <typename T>
    const T& getValueOrDefault(MDLabel label, const T& def) const { return MDRow::getValueOrDefault<T>(label, def); }

    template <typename T>
    T& getValueOrDefault(MDLabel label, const T& def) { return MDRow::getValueOrDefault<T>(label, def); }

    template <typename T, typename T1>
    void getValueOrDefault(MDLabel label, T &d, T1 def) const { // FIXME: deprecated
        return MDRow::getValueOrDefault<T, T1>(label, d, def);
    }

    template <typename T>
    void setValue(MDLabel label, const T &d, bool addLabel = true) { return MDRow::setValue<T>(label, d, addLabel); }
};

#endif
