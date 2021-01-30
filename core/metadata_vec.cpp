/**************************************************************************
 *
 * Authors:      J.R. Bilbao-Castro (jrbcast@ace.ual.es)
 *               Jan Horacek (xhorace4@fi.muni.cz)
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

#include "metadata_vec.h"

MetaDataVec::MetaDataVec() {
    init(std::vector<MDLabel>());
}

MetaDataVec::MetaDataVec(const std::vector<MDLabel> &labelsVector) {
    init(labelsVector);
}

MetaDataVec::MetaDataVec(const FileName &fileName, const std::vector<MDLabel> &desiredLabels) {
    init(desiredLabels);
    read(fileName);
}

void MetaDataVec::init(const std::vector<MDLabel> &labelsVector) {
    size_t col = 0;
    for (const auto label : labelsVector) {
        _label_to_col[label] = col;
        col++;
    }
    _no_columns = col;
}

void MetaDataVec::read(const FileName &inFile) {
    throw NotImplemented();
}

void MetaDataVec::write(const FileName &outFile, WriteModeMetaData mode) const {
    String blockName;
    FileName extFile;
    FileName _outFile;

    blockName = outFile.getBlockName();
    if (blockName.empty())
        blockName = DEFAULT_BLOCK_NAME;
    _outFile = outFile.removeBlockName();
    extFile = outFile.getExtension();

    if (extFile == "xml") {
        writeXML(_outFile, blockName, mode);
    } else if (extFile == "sqlite") {
        throw NotImplemented();
    } else {
        writeStar(_outFile, blockName, mode);
    }
}

size_t MetaDataVec::addRow(const MDRowVec &row) {
    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(label);
        size_t rowCol = (*row._label_to_col)[label];
        if (rowCol > -1)
            if (!this->containsLabel(label))
                this->addLabel(label);
    }

    MetaDataVecRow newRow;

    for (size_t coli = 0; coli < _no_columns; coli++)
        newRow.push_back(_col_to_label[coli]);

    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(label);
        size_t rowCol = (*row._label_to_col)[label];
        if (rowCol > -1) {
            size_t ourCol = _label_to_col[label];
            newRow[ourCol] = (*row._row)[rowCol];
        }
    }

    _rows.push_back(newRow);
}
