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

/** Class for holding an entire row of posible MDObject */
/* Row could be attached to OR detached from metadata
 * - Detached row: holds its own _row & _label_to_col. When copied, content of
 *   these objects is copied, not just pointer.
 * - Attached row: _row and _label_to_col points directly to Metadata object.
 *   When copied, just pointers are copied. Assumed to be read-only.
 */
class MDRowVec : public MDRow {
private:
    std::vector<MDObject>* _row;
    size_t _rowi;
    std::array<int, MDL_LAST_LABEL>* _label_to_col;
    bool _in_metadata;

    MDObject* iteratorValue(size_t i) const override;

public:
    MDRowVec();
    MDRowVec(std::vector<MDObject>& row, size_t rowi, std::array<int, MDL_LAST_LABEL>& label_to_col);
    MDRowVec(const MDRowVec&);
    MDRowVec& operator = (const MDRowVec&);
    virtual ~MDRowVec();

    bool empty() const override;
    int size() const override;
    void clear() override;
    bool inMetadata() const;

    bool containsLabel(MDLabel label) const override;
    std::vector<MDLabel> labels() const override;
    void addLabel(MDLabel label) override;

    MDObject *getObject(MDLabel label) const override;

    bool getValue(MDObject &object) const override;
    void setValue(const MDObject &object) override;

    friend std::ostream& operator << (std::ostream &out, const MDRowVec &row);

    // Templated functions from based class must be retemplated

    template <typename T>
    bool getValue(MDLabel label, T &d) const { return MDRow::getValue(label, d); }

    template <typename T, typename T1>
    void getValueOrDefault(MDLabel label, T &d, T1 def) const { return MDRow::getValueOrDefault(label, d, def); }

    template <typename T>
    void setValue(MDLabel label, const T &d, bool addLabel = true) { return MDRow::setValue(label, d, addLabel); }

    friend class MetaDataVec;
};

#endif
