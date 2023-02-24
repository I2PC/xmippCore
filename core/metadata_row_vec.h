/***************************************************************************
 *
 * Authors:     Jan Horacek (xhorace4@fi.muni.cz)
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

#ifndef CORE_METADATA_ROW_VEC_H
#define CORE_METADATA_ROW_VEC_H

#include "metadata_row_base.h"

#include <array>

/** Class for holding an entire row of MDObject in MetaDataVec.
 * Row could be attached to OR detached from metadata
 * - Detached row: holds its own _row & _label_to_col. When copied, content of
 *   these objects is copied, not just pointer.
 * - Attached row: _row and _label_to_col points directly to Metadata object.
 *   When copied, just pointers are copied. Assumed to be read-only.
 *
 * ### Notes
 *  1. It's fast to create MDRowVec from MetaDataVec row, because only pointers
 *     are initialized. No data are copied.
 *  2. Adding MDRowVec to MetaDataVec requires copying whole row.
 */
class MDRowVec : public MDRow {
private:
    std::vector<MDObject>* _row;
    size_t _rowi;
    std::array<int, MDL_LAST_LABEL>* _label_to_col;
    std::vector<MDLabel>* _col_to_label;
    size_t* _no_columns; // global number of columns in whole MetaData
    bool _in_metadata;

    MDObject* iteratorValue(size_t i) override;
    const MDObject* iteratorValue(size_t i) const override;
    size_t newCol(const MDLabel);

public:
    static MDRowVec deepCopy(const MDRowVec&);

    MDRowVec();
    MDRowVec(std::vector<MDObject>& row, size_t rowi, std::array<int, MDL_LAST_LABEL>& label_to_col,
             std::vector<MDLabel>& col_to_label, size_t& no_columns);
    MDRowVec(const std::vector<MDObject>& row, size_t rowi, const std::array<int, MDL_LAST_LABEL>& label_to_col,
             const std::vector<MDLabel>& col_to_label, const size_t& no_columns);
    MDRowVec(const MDRowVec&);
    MDRowVec& operator = (const MDRowVec&);
    MDRow& operator = (const MDRow& row);
    virtual ~MDRowVec();

    void detach() override;

    bool empty() const override;
    int size() const override;
    void clear() override;
    bool inMetadata() const;

    bool containsLabel(MDLabel label) const override;
    std::vector<MDLabel> labels() const override;
    void addLabel(MDLabel label) override;

    MDObject *getObject(MDLabel label) override;
    const MDObject *getObject(MDLabel label) const override;

    void setValue(const MDObject &object) override;

    friend std::ostream& operator << (std::ostream &out, const MDRowVec &row);

    // Templated functions from based class must be retemplated

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

    friend class MetaDataVec;
};

#endif
