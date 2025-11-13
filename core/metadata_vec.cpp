/**************************************************************************
 i
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

#include <algorithm>
#include <cassert>
#include <fstream>
#include <random>
#include "metadata_vec.h"
#include "metadata_generator.h"
#include "xmipp_image.h"

#include <sys/stat.h>
#include <fcntl.h>
#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif
#include "xmipp_funcs.h"

MetaDataVec::MetaDataVec() {
    init({});
}

MetaDataVec::MetaDataVec(const std::vector<MDLabel> &labelsVector) {
    init(labelsVector);
}

MetaDataVec::MetaDataVec(const FileName &fileName, const std::vector<MDLabel> &desiredLabels) {
    init(desiredLabels);
    read(fileName);
}

MetaDataVec::MetaDataVec(const FileName &fileName) {
    init({});
    read(fileName);
}

MetaDataVec::MetaDataVec(const MetaData &md) {
    init({});
    MetaData::operator=(md);
}


void MetaDataVec::init(const std::vector<MDLabel> &labelsVector) {
    this->clear();
    _label_to_col[MDL_OBJID] = 0;
    _col_to_label.push_back(MDL_OBJID);
    size_t col = 1;
    for (const auto label : labelsVector) {
        if (label != MDL_OBJID) {
            _label_to_col[label] = col;
            _col_to_label.push_back(label);
            col++;
        }
    }
    _no_columns = col;
}

int MetaDataVec::_labelIndex(MDLabel label) const {
    return this->_label_to_col[label];
}

const MDObject& MetaDataVec::_getObject(size_t i, MDLabel label) const {
    return this->_getObject(this->_rows.at(i), label);
}

MDObject& MetaDataVec::_getObject(size_t i, MDLabel label) {
    return this->_getObject(this->_rows.at(i), label);
}

const MDObject& MetaDataVec::_getObject(const MetaDataVecRow& row, MDLabel label) const {
    int labelIndex = this->_labelIndex(label);
    if ((labelIndex < 0) || (static_cast<size_t>(labelIndex) >= row.size()))
        throw ColumnDoesNotExist(label, getFilename());
    return row.at(labelIndex);
}

MDObject& MetaDataVec::_getObject(MetaDataVecRow& row, MDLabel label) const {
    int labelIndex = this->_labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist(label, getFilename());
    return row.at(labelIndex);
}

int MetaDataVec::_rowIndex(size_t id) const {
    if (this->_id_to_index.find(id) == this->_id_to_index.end())
        return -1;
    return this->_id_to_index.at(id);
}

size_t MetaDataVec::_rowIndexSafe(size_t id) const {
    int i = this->_rowIndex(id);
    if (i == -1)
        throw ObjectDoesNotExist(id, getFilename());
    return i;
}

void MetaDataVec::read(const FileName &filename, const std::vector<MDLabel> *desiredLabels, bool decomposeStack) {
    String blockName;
    FileName inFile;

    blockName = filename.getBlockName();
    inFile = filename.removeBlockName();
    String extFile = filename.getExtension();
    blockName = escapeForRegularExpressions(blockName);

    this->clear();
    this->setColumnFormat(true);

    if (extFile == "xml")
        this->readXML(_inFile, desiredLabels, blockName, decomposeStack);
    else if (extFile == "sqlite")
        throw NotImplemented("Reading from .sqlite file into MetaDataVec not implemented!");
    else
        this->readStar(filename, desiredLabels, blockName, decomposeStack);
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
        throw NotImplemented("Writing to .sqlite file from MetaDataVec not implemented!");
    } else {
        writeStar(_outFile, blockName, mode);
    }
}

void MetaDataVec::writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const {
    // FIXME: implement
    throw NotImplemented("writeXML not implemented");
}

void MetaDataVec::writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const {
    // FIXME: implement
    throw NotImplemented("writeText not immplemented");
}

void MetaDataVec::clear() {
    MetaData::clear();
    this->_rows.clear();
    std::fill(this->_label_to_col.begin(), this->_label_to_col.end(), -1);
    this->_col_to_label.clear();
    this->_id_to_index.clear();
    this->_no_columns = 0;
    this->_next_id = 1;
}

void MetaDataVec::_setRow(const MDRow &row, size_t index) {
    if (dynamic_cast<const MDRowVec*>(&row) != nullptr) {
        // No not change same row
        const MDRowVec& mdRowVec = dynamic_cast<const MDRowVec&>(row);
        if ((mdRowVec._in_metadata) && (mdRowVec._rowi == index) && (mdRowVec._row == &this->_rows[index]))
            return;
    }

    size_t newRowSize = 0;
    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(labeli);
        if (row.containsLabel(label)) {
            if (!this->containsLabel(label))
                newRowSize = this->size();
            else if (this->_label_to_col[label]+1 > newRowSize)
                newRowSize = this->_label_to_col[label]+1;
        }
    }

    for (size_t column = 0; (column < this->_col_to_label.size()) && (column < newRowSize); ++column)
        if ((this->_col_to_label[column] != MDL_OBJID) && (!row.containsLabel(this->_col_to_label[column])))
            throw ColumnDoesNotExist("New row does not contain required MetaData column: "+
                                     MDL::label2Str(this->_col_to_label[column])+"!");

    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(labeli);
        if (row.containsLabel(label) && !this->containsLabel(label))
            this->addLabel(label);
    }

    MetaDataVecRow& editRow = this->_rows[index];
    editRow.clear();

    for (size_t coli = 0; coli < this->_no_columns; coli++)
        editRow.push_back({this->_col_to_label[coli]});

    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(labeli);
        if (row.containsLabel(label)) {
            size_t ourCol = _label_to_col[label];
            row.getValue(editRow[ourCol]);
        }
    }
}

size_t MetaDataVec::addRow(const MDRow &row) {
    /* Id:
     * When ‹row› does not contain MDL_OBJID column, it is created from ‹this->_nextId›
     * When ‹row› contains id which is NOT present in Metadata, id is kept.
     * When ‹row› contains id which IS present in Metadata, id is changed from ‹this->_nextId›.
     * When ‹this->_nextId› is present in MetaData, assert fails. This should not happen.
    */

    MetaDataVecRow newRow;
    _rows.emplace_back(newRow);
    this->_setRow(row, _rows.size()-1);

    if ((!row.containsLabel(MDL_OBJID)) || (row.getValue<size_t>(MDL_OBJID) == BAD_OBJID)) {
        MetaDataVecRow& _row = this->_rows[_rows.size()-1];
        if (!this->containsLabel(MDL_OBJID))
            this->addLabel(MDL_OBJID);
        this->_expand(_row, MDL_OBJID);
        _row[this->_labelIndex(MDL_OBJID)] = MDObject(MDL_OBJID, this->_next_id);
    }

    size_t rowId = getRowId(_rows.size()-1);

    if (this->_id_to_index.find(rowId) != this->_id_to_index.end()) {
        MetaDataVecRow& _row = this->_rows[_rows.size()-1];
        _row[this->_labelIndex(MDL_OBJID)] = MDObject(MDL_OBJID, this->_next_id);
        rowId = this->_next_id;
    }

    assert(this->_id_to_index.find(rowId) == this->_id_to_index.end());

    if (rowId >= this->_next_id)
        this->_next_id = rowId+1;
    this->_id_to_index[rowId] = _rows.size()-1;
    return rowId;
}

void MetaDataVec::addRows(const std::vector<MDRowVec> &rows) {
    for (const auto row : rows)
        this->addRow(row);
}

int MetaDataVec::getMaxStringLength(const MDLabel thisLabel) const {
    return 255;
}

bool MetaDataVec::setValueCol(const MDObject &mdValueIn) {
    const auto &label = mdValueIn.label;
    int labelIndex = this->_labelIndex(label);
    if (labelIndex < 0) {
        this->addLabel(label);
        labelIndex = this->_labelIndex(label);
        for (auto &r : _rows) _expand(r, labelIndex);
    }
    for (auto &r : _rows) r[labelIndex] = mdValueIn;
    return true;
}

bool MetaDataVec::setValue(const MDObject &mdValueIn, size_t id) {
    if (!this->containsLabel(mdValueIn.label))
        this->addLabel(mdValueIn.label);

    MetaDataVecRow& row = this->_rows[this->_rowIndexSafe(id)];
    this->_expand(row, mdValueIn.label);

    row[this->_labelIndex(mdValueIn.label)] = mdValueIn;
    return true;
}

bool MetaDataVec::getValue(MDObject &mdValueOut, size_t id) const {
    try {
        mdValueOut = this->_getObject(this->_rowIndexSafe(id), mdValueOut.label);
    } catch (const ColumnDoesNotExist&) {
        return false;
    }
    return true;
}

MDObject &MetaDataVec::getValue(MDLabel label, size_t id) {
    return this->_getObject(this->_rowIndexSafe(id), label);
}

const MDObject &MetaDataVec::getValue(MDLabel label, size_t id) const {
    return this->_getObject(this->_rowIndexSafe(id), label);
}

std::unique_ptr<MDRow> MetaDataVec::getRow(size_t id) {
    int i = this->_rowIndex(id);
    if (i < 0)
        return nullptr;
    return memoryUtils::make_unique<MDRowVec>(
        this->_rows[i], i, this->_label_to_col, this->_col_to_label, this->_no_columns
    );
}

std::unique_ptr<const MDRow> MetaDataVec::getRow(size_t id) const {
    int i = this->_rowIndex(id);
    if (i < 0)
        return nullptr;
    return memoryUtils::make_unique<MDRowVec>(
        this->_rows[i], i, this->_label_to_col, this->_col_to_label, this->_no_columns
    );
}

MDRowVec MetaDataVec::getRowVec(size_t id) {
    size_t i = this->_rowIndexSafe(id);
    return MDRowVec(this->_rows[i], i, this->_label_to_col, this->_col_to_label, this->_no_columns);
}

const MDRowVec MetaDataVec::getRowVec(size_t id) const {
    size_t i = this->_rowIndexSafe(id);
    return MDRowVec(this->_rows[i], i, this->_label_to_col, this->_col_to_label, this->_no_columns);
}

void MetaDataVec::getRow(MDRowVec &row, size_t id) {
    size_t i = this->_rowIndexSafe(id);
    row = MDRowVec(this->_rows[i], i, this->_label_to_col, this->_col_to_label, this->_no_columns);
}

bool MetaDataVec::getRowValues(size_t id, std::vector<MDObject> &values) const {
    int i = this->_rowIndex(id);
    if (i < 0)
        return false;
    values = this->_rows[i];
    return true;
}

size_t MetaDataVec::getRowId(size_t i) const {
    return this->_getObject(i, MDL_OBJID).getValue2(size_t());
}

size_t MetaDataVec::getRowId(const MetaDataVecRow& row) const {
    int labelIndex = _labelIndex(MDL_OBJID);
    if (labelIndex < 0)
        throw ColumnDoesNotExist(MDL_OBJID, getFilename());
    return row.at(labelIndex).getValue2(size_t());
}

void MetaDataVec::getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const {
    valuesOut.clear();
    int labelIndex = this->_labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist(label, getFilename());
    for (const auto& vecRow : this->_rows)
        valuesOut.emplace_back(vecRow.at(labelIndex));
}

void MetaDataVec::setColumnValues(const std::vector<MDObject> &valuesIn) {
    for (size_t i = 0; i < std::min(valuesIn.size(), this->_rows.size()); i++) {
        int labelIndex = this->_labelIndex(valuesIn[i].label);
        if (labelIndex < 0)
            this->addLabel(valuesIn[i].label);
        labelIndex = this->_labelIndex(valuesIn[i].label);
        if (labelIndex < 0)
            throw ColumnDoesNotExist(valuesIn[i].label, getFilename());
        this->_expand(this->_rows[i], valuesIn[i].label);
        this->_rows[i][labelIndex] = valuesIn[i];
    }
}

bool MetaDataVec::setRow(const MDRow &row, size_t id) {
    this->_setRow(row, this->_rowIndexSafe(id));
    return true;
}

bool MetaDataVec::isEmpty() const {
    return this->_rows.empty();
}

size_t MetaDataVec::size() const {
    return this->_rows.size();
}

bool MetaDataVec::containsLabel(const MDLabel label) const {
    return this->_labelIndex(label) > -1;
}

bool MetaDataVec::addLabel(const MDLabel label, int pos) {
    if (pos != -1)
        throw NotImplemented("addLabel to -1 not implemented");
    if (this->_label_to_col[label] != -1)
        return true;

    this->_no_columns++;
    size_t column = this->_no_columns-1;
    this->_label_to_col[label] = column;
    this->_col_to_label.emplace_back(label);
    return true;
}

bool MetaDataVec::removeLabel(const MDLabel label) {
    if (this->_label_to_col[label] == -1)
        return false;
    int column = this->_label_to_col[label];

    for (auto& vecRow : this->_rows)
        if (static_cast<int>(vecRow.size()) > column)
            vecRow.erase(vecRow.begin()+column); // this is expensive

    // FIXME: test this properly
    this->_label_to_col[label] = -1;
    for (size_t i = 0; i < MDL_LAST_LABEL; i++) {
        if (this->_label_to_col[i] > column) {
            this->_label_to_col[i]--;
            this->_col_to_label[this->_label_to_col[i]] = MDLabel(i);
        }
    }
    this->_no_columns--;
    return true;
}

size_t MetaDataVec::addObject() {
    MDRowVec row;
    row.setValue(MDObject(MDL_OBJID, this->_next_id));
    return this->addRow(row);
}

void MetaDataVec::importObject(const MetaData &md, const size_t id, bool doClear) {
    if (doClear) {
        this->clear();
        this->copyInfo(md);
    }

    std::unique_ptr<const MDRow> row = md.getRow(id);
    if (row == nullptr)
        return;

    MDRowVec newRow;
    for (const MDObject* obj : *row)
        if (obj->label != MDL_OBJID)
            newRow.setValue(*obj);
    this->addRow(newRow);
}

void MetaDataVec::importObjects(const MetaData &md, const std::vector<size_t> &objectsToAdd, bool doClear) {
    if (doClear) {
        this->clear();
        this->copyInfo(md);
    }

    for (size_t objId : objectsToAdd)
        this->importObject(md, objId, false);
}

void MetaDataVec::importObjects(const MetaData &md, const MDQuery &query, bool doClear) {
    // FIXME: move this to MetaData?
    std::vector<size_t> ids;
    md.findObjects(ids, query);
    this->importObjects(md, ids, doClear);
}

bool MetaDataVec::removeObject(size_t id) {
    int i = this->_rowIndex(id);
    if (i < 0)
        return false;

    this->_id_to_index.erase(id);
    this->_rows.erase(this->_rows.begin()+i);

    for (size_t j = i; j < this->_rows.size(); j++)
        this->_id_to_index[this->getRowId(j)]--;

    return true;
}

void MetaDataVec::removeObjects(const std::vector<size_t> &toRemove) {
    for (size_t id : toRemove)
        this->removeObject(id);
}

int MetaDataVec::removeObjects() {
    size_t count = this->size();
    this->clear();
    return count;
}

int MetaDataVec::removeObjects(const MDQuery& query) {
    // FIXME: move this to MetaData?
    std::vector<size_t> ids;
    this->findObjects(ids, query);
    this->removeObjects(ids);
    return true;
}

size_t MetaDataVec::firstRowId() const {
    return this->getRowId(0);
}

size_t MetaDataVec::lastRowId() const {
    return this->getRowId(this->size()-1);
}

bool MetaDataVec::_match(const MetaDataVecRow& row, const MDQuery& query) const {
    if (dynamic_cast<const MDMultiQuery*>(&query) != nullptr) {
        // Process MDMultiQuery
        const MDMultiQuery& mq = dynamic_cast<const MDMultiQuery&>(query);
        if (mq.operations.size() == 0)
            return true; // emptuy query matches always
        String operation = mq.operations[0];
        for (const auto& op : mq.operations)
            assert(op == operation); // support only all same operations
        for (const auto& q : mq.queries) {
            if ((operation == "AND") && (!this->_match(row, *q)))
                return false;
            else if ((operation == "OR") && (this->_match(row, *q)))
                return true;
        }
        return (operation == "AND"); // false for "OR", true for "AND"
    }

    if (dynamic_cast<const MDValueRelational*>(&query) == nullptr)
        throw NotImplemented("_match for this type of query not implemented");
            // MDValueRange, MDExpression, MDMultiQuery not implemented yet
            // MDExpression will probably never be supported as it is raw SQL expression

    const MDValueRelational& rel = dynamic_cast<const MDValueRelational&>(query);

    if (rel.value == nullptr)
        return false;

    size_t labeli = this->_labelIndex(rel.value->label);
    if (labeli >= row.size())
        return false;

    const MDObject& mdObj = row[labeli];

    if (rel.op == RelationalOp::EQ)
        return *(rel.value) == mdObj;
    if (rel.op == RelationalOp::NE)
        return *(rel.value) != mdObj;
    if (rel.op == RelationalOp::GT) // FIXME: check if < & > are not swapped
        return *(rel.value) < mdObj;
    if (rel.op == RelationalOp::GE)
        return *(rel.value) <= mdObj;
    if (rel.op == RelationalOp::LT)
        return *(rel.value) > mdObj;
    if (rel.op == RelationalOp::LE)
        return *(rel.value) >= mdObj;

    throw std::logic_error("MetaDataVec::_match: unknown operator");
}

size_t MetaDataVec::firstObject(const MDQuery& query) const {
    // FIXME: should first be first in _rows order or ids order?
    for (const MetaDataVecRow& row : this->_rows)
        if (this->_match(row, query))
            return this->getRowId(row);
    return BAD_OBJID;
}

void MetaDataVec::findObjects(std::vector<size_t> &objectsOut, const MDQuery &query) const {
    // FIXME: should first be first in _rows order or ids order?
    objectsOut.clear();
    for (const MetaDataVecRow& row : this->_rows)
        if (this->_match(row, query))
            objectsOut.emplace_back(this->getRowId(row));
}

void MetaDataVec::findObjects(std::vector<size_t> &objectsOut, int limit) const {
    objectsOut.clear();
    for (size_t i = 0; i < std::min<size_t>(limit, this->size()); i++)
        objectsOut.emplace_back(this->getRowId(this->_rows[i]));
}

size_t MetaDataVec::countObjects(const MDQuery& query) const {
    size_t count = 0;
    for (const MetaDataVecRow& row : this->_rows)
        if (this->_match(row, query))
            count++;
    return count;
}

bool MetaDataVec::containsObject(size_t objectId) const {
    return this->_rowIndex(objectId) > -1;
}

bool MetaDataVec::containsObject(const MDQuery& query) const {
    for (const MetaDataVecRow& row : this->_rows)
        if (this->_match(row, query))
            return true;
    return false;
}

bool MetaDataVec::containsObject(size_t objectId) {
    return this->_id_to_index.find(objectId) != this->_id_to_index.end();
}


void MetaDataVec::_writeRows(std::ostream &os) const {
    for (const MetaDataVecRow& row : this->_rows) {
        for (size_t i = 0; i < MDL_LAST_LABEL; i++) {
            const MDLabel label = static_cast<MDLabel>(i);
            if ((label != MDL_STAR_COMMENT) && (label != MDL_OBJID) && (this->_label_to_col[i] > -1)) {
                os.width(1);
                if (this->_labelIndex(label) < static_cast<int>(row.size()))
                    this->_getObject(row, label).toStream(os, true);
                else
                    throw ColumnDoesNotExist(label, getFilename());
                os << " ";
            }
        }
        os << '\n';
    }
}

void MetaDataVec::write(std::ostream &os, const String &blockName, WriteModeMetaData mode) const {
    if (mode==MD_OVERWRITE)
        os << FileNameVersion << " * "// << (isColumnFormat ? "column" : "row")
        << '\n' // write which type of format (column or row) and the path;
        << WordWrap(this->_comment, line_max);     // write md comment in the 2nd comment line of header

    // write data block
    String _szBlockName("data_");
    _szBlockName += blockName;

    if (this->isColumnFormat()) {
        // write md columns in 3rd comment line of the header
        os << _szBlockName << '\n';
        os << "loop_" << '\n';
        for (size_t i = 0; i < MDL_LAST_LABEL; i++)
            if ((i != MDL_STAR_COMMENT) && (i != MDL_OBJID) && (this->_label_to_col[i] > -1))
                os << " _" << MDL::label2Str(static_cast<MDLabel>(i)) << '\n';
        _writeRows(os);

        //Put the activeObject to the first, if exists
    } else { // row format
        os << _szBlockName << '\n';

        // Print single object
        assert(this->_rows.size() == 1);

        for (size_t i = 0; i < MDL_LAST_LABEL; i++) {
            const MDLabel label = static_cast<MDLabel>(i);
            if ((label != MDL_STAR_COMMENT) && (label != MDL_OBJID) && (this->_label_to_col[i] > -1)) {
                os << " _" << MDL::label2Str(label) << " ";
                if (this->_labelIndex(label) < static_cast<int>(this->_rows[0].size()))
                    this->_getObject(this->_rows[0], label).toStream(os);
                else
                    throw ColumnDoesNotExist(label, getFilename());
                os << '\n';
            }
        }
    }
}

void MetaDataVec::_parseObjects(std::istream &is, std::vector<MDObject*> &columnValues,
                                const std::vector<MDLabel> *desiredLabels, bool firstTime) {
    for (size_t i = 0; i < columnValues.size(); i++) {
        columnValues[i]->fromStream(is);

        if (is.fail()) {
           String errorMsg = formatString("MetaData: Error parsing column '%s' value.",
                                          MDL::label2Str(columnValues[i]->label).c_str());
           columnValues[i]->failed = true;
           std::cerr << "WARNING: " << errorMsg << std::endl;
           //REPORT_ERROR(ERR_MD_BADLABEL, (String)"read: Error parsing data column, expecting " + MDL::label2Str(object.label));
        } else if (firstTime) {
            // Check if current column label exists.
            if (columnValues[i]->label != MDL_UNDEFINED) {
                // If there are no desired labels then add all.
                bool reallyAdd=false;
                if (desiredLabels == NULL) {
                    reallyAdd = true;
                } else {
                    // Check if current column belongs to desired labels.
                    for (size_t j = 0; j < desiredLabels->size(); ++j) {
                        if ((*desiredLabels)[j] == columnValues[i]->label) {
                            reallyAdd = true;
                            break;
                        }
                    }
                }

                // Add label if not exists.
                if (reallyAdd)
                    this->addLabel(columnValues[i]->label);
            }
        }
    }

    // Insert elements in DB.
    MDRowVec newRow;
    for (size_t i = 0; i < columnValues.size(); i++)
        if (columnValues[i] != nullptr)
            newRow.setValue(*columnValues[i]);
    this->addRow(newRow);
}

void MetaDataVec::readXML(const FileName &inFile, const std::vector<MDLabel> *desiredLabels, const String & blockRegExp, bool decomposeStack) {
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "readXML not implemented yet");
}

void MetaDataVec::readPlain(const FileName &inFile, const String &labelsString, const String &separator) {
    // TODO
    throw NotImplemented("readPlain not implemented");
}

void MetaDataVec::addPlain(const FileName &inFile, const String &labelsString, const String &separator) {
    // TODO
    throw NotImplemented("addPlain not implemented");
}

double MetaDataVec::getColumnMax(MDLabel column) {
    // TODO
    throw NotImplemented("getColumnMax not implemented");
}

double MetaDataVec::getColumnMin(MDLabel column) {
    // TODO
    throw NotImplemented("getColumnMin not implemented");
}

void MetaDataVec::replace(const MDLabel label, const String &oldStr, const String &newStr) {
    // TODO
    throw NotImplemented("replace not implemented");
}

void MetaDataVec::randomize(const MetaData &MDin) {
    *this = MDin;
    std::random_device rd;
    auto g = std::mt19937(rd());
    std::shuffle(this->_rows.begin(), this->_rows.end(), g);
    this->_recalc_id_to_index();
}

void MetaDataVec::_recalc_id_to_index() {
    this->_id_to_index.clear();
    for (size_t i = 0; i < this->_rows.size(); i++)
        this->_id_to_index[this->getRowId(i)] = i;
}

void MetaDataVec::removeDuplicates(MetaData &MDin, MDLabel label) {
    *this = MDin; // FIXME: maybe join?
    std::vector<MetaDataVecRow> new_rows;
    for (const auto& row : this->_rows) {
        if (!this->_contains(new_rows, row))
            new_rows.emplace_back(row);
    }
    this->_rows = new_rows;
}

bool MetaDataVec::_contains(const std::vector<MetaDataVecRow>& rows, const MetaDataVecRow& row) const {
    for (const auto& _row : rows)
        if (this->_rowsEq(row, _row))
            return true;
    return false;
}

bool MetaDataVec::_rowsEq(const MetaDataVecRow& a, const MetaDataVecRow& b) const {
    for (size_t label = 0; label < MDL_LAST_LABEL; label++) {
        if ((label == MDL_COMMENT) || (label == MDL_OBJID))
            continue;

        int labeli = this->_labelIndex(static_cast<MDLabel>(label));
        if (labeli > -1) { // label is active
            if ((static_cast<size_t>(labeli) < a.size()) !=
                (static_cast<size_t>(labeli) < b.size()))
                return false; // item present in one row, but not other
            if (static_cast<size_t>(labeli) >= a.size())
                continue; // label not present in both rows
            if (!a[labeli].eq(b[labeli], this->precision()))
                return false; // MDObjects are diffrent
        }
    }
    return true;
}

void MetaDataVec::sort(const MetaDataVec &MDin, const MDLabel sortLabel, bool asc, int limit, int offset) {
    *this = MDin;

    int label_index = this->_labelIndex(sortLabel);
    if (label_index > -1) {
        std::sort(this->_rows.begin(), this->_rows.end(),
            [label_index, asc](const MetaDataVecRow &a, const MetaDataVecRow &b) {
                if (asc)
                    return a[label_index] < b[label_index];
                return a[label_index] > b[label_index];
            }
        );

        this->_rows.erase(this->_rows.begin(), this->_rows.begin()+offset);
        if ((limit > 0) && (limit < this->_rows.size()))
            this->_rows.erase(this->_rows.begin()+limit, this->_rows.end());

        this->_recalc_id_to_index();

    }
}

void MetaDataVec::sort(MetaDataVec &MDin, const String &sortLabel, bool asc, int limit, int offset) {
    // TODO
    throw NotImplemented("sort not implemented");
}

void MetaDataVec::split(size_t parts, std::vector<MetaDataVec> &results, const MDLabel sortLabel) const {
    if (parts > this->size())
        REPORT_ERROR(ERR_MD, "MetaDataDb::split: Couldn't split a metadata in more parts than its size");

    MetaDataVec sorted;
    if (sortLabel == MDL_UNDEFINED)
        sorted = *this;
    else
        sorted.sort(*this, sortLabel);

    results.clear();
    results.resize(parts);
    for (size_t i = 0; i < parts; i++) {
        MetaDataVec &md = results[i];
        size_t firsti, lasti;
        divide_equally(sorted.size(), parts, i, firsti, lasti);
        for (size_t j = firsti; j <= lasti; j++)
            md.addRow(sorted.getRowVec(sorted.getRowId(j)));
    }
}

void MetaDataVec::selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel) {
    // TODO
    throw NotImplemented("selectRandomSubset not implemented");
}

void MetaDataVec::selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                             const MDLabel sortLabel) {
    this->sort(mdIn, sortLabel, true, numberOfObjects, startPosition);
}

/*void makeAbsPath(const MDLabel label=MDL_IMAGE);*/


void MetaDataVec::fillConstant(MDLabel label, const String &value) {
    // FIXME: move to MetaData and use common MDGenerator?
    MDConstGenerator generator(value);
    generator.label = label;
    generator.fill(*this);
}

void MetaDataVec::fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3) {
    // FIXME: move to MetaData and use common MDGenerator?
    MDRandGenerator generator(op1, op2, mode, op3);
    generator.label = label;
    generator.fill(*this);
}

void MetaDataVec::fillLinear(MDLabel label, double initial, double step) {
    // FIXME: move to MetaData and use common MDGenerator?
    MDLinealGenerator generator(initial, step);
    generator.label = label;
    generator.fill(*this);
}

void MetaDataVec::_expand(MetaDataVecRow& row, const MDLabel label) {
    int labeli = this->_labelIndex(label);
    if (labeli < 0)
        this->addLabel(label);
    this->_expand(row,  this->_labelIndex(label));
}

void MetaDataVec::_expand(MetaDataVecRow& row, size_t labeli) {
    // In assert: all labels to labeli (including) must be present in
    // this->_col_to_label.

    if (labeli < row.size())
        return; // space for label already present

    for (size_t i = row.size(); i <= labeli; i++)
        row.emplace_back(MDObject(this->_col_to_label.at(i)));
}

void MetaDataVec::copyColumn(MDLabel labelDest, MDLabel labelSrc) {
    int labelsrci = this->_labelIndex(labelSrc);
    int labeldesti = this->_labelIndex(labelDest);
    if (labelsrci < 0)
        return;
    if (labeldesti < 0)
        this->addLabel(labelDest);
    labeldesti = this->_labelIndex(labelDest);

    for (MetaDataVecRow& row : this->_rows) {
        if (static_cast<size_t>(labeldesti) >= row.size())
            this->_expand(row, labeldesti);

        if (static_cast<size_t>(labelsrci) < row.size()) {
            row[labeldesti] = row[labelsrci];
            row[labeldesti].label = labelDest;
        }
        // else row[labeldesti] is empty MDObject (from previous if)
    }
}

void MetaDataVec::copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc) {
    // TODO
    throw NotImplemented("copyColumnTo not implemented");
}

void MetaDataVec::renameColumn(MDLabel oldLabel, MDLabel newLabel) {
    assert(!this->containsLabel(newLabel));
    int labeloldi = this->_labelIndex(oldLabel);
    if (labeloldi < 0)
        throw ColumnDoesNotExist(oldLabel, getFilename());

    this->_label_to_col[newLabel] = labeloldi;
    this->_label_to_col[oldLabel] = -1;
    this->_col_to_label[labeloldi] = newLabel;

    for (auto& row : this->_rows) {
        if (labeloldi < static_cast<int>(row.size()))
            row[labeloldi].label = newLabel;
    }
}

void MetaDataVec::renameColumn(const std::vector<MDLabel> &oldLabel,
                               const std::vector<MDLabel> &newLabel) {
    // TODO
    throw NotImplemented("renameColumn not implemented");
}


bool MetaDataVec::operator==(const MetaDataVec& op) const {
    // This comparison ignores order of labels and row ids, everything else must be same.

    if (this->_rows.size() != op._rows.size())
        return false;

    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; labeli++)
        if ((this->_label_to_col[labeli] > -1) != (op._label_to_col[labeli] > -1))
            return false;

    for (size_t i = 0; i < this->_rows.size(); i++) {
        for (size_t labeli = 0; labeli < MDL_LAST_LABEL; labeli++) {
            if ((labeli == MDL_COMMENT) || (labeli == MDL_OBJID))
                continue;

            int thisLabelColI = this->_label_to_col[labeli];
            int opLabelColI = op._label_to_col[labeli];
            if (thisLabelColI > -1) {
                if ((static_cast<size_t>(thisLabelColI) < this->_rows[i].size()) !=
                    (static_cast<size_t>(opLabelColI) < op._rows[i].size()))
                    return false; // item present in one row, but not other
                if (static_cast<size_t>(thisLabelColI) >= this->_rows[i].size())
                    continue; // label not present in both rows
                if (!this->_rows[i][thisLabelColI].eq(op._rows[i][opLabelColI], this->precision()))
                    return false; // MDObjects are diffrent
            }
        }
    }

    return true; // all rows same → ok
}

std::vector<MDLabel> MetaDataVec::getActiveLabels() const {
    std::vector<MDLabel> out;
    for (size_t i = MDL_GATHER_ID; i < MDL_LAST_LABEL; i++) // ignore MDL_FIRST_LABEL = MDL_OBJID
        if (this->_label_to_col[i] > -1)
            out.emplace_back(static_cast<MDLabel>(i));
    return out;
}

std::ostream& operator<<(std::ostream& o, const MetaDataVec& md) {
    md.write(o);
    return o;
}
