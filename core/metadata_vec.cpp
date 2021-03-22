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

#include <algorithm>
#include "metadata_vec.h"
#include "metadata_generator.h"

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

MetaDataVec::MetaDataVec(const FileName &fileName) {
    init({});
    read(fileName);
}

MetaDataVec::MetaDataVec(const MetaData &md) {
    *this = md;
}

MetaDataVec& MetaDataVec::operator=(const MetaData &md) {
    this->copyInfo(md);
    for (const auto& row : md)
        this->addRow(row);
    return *this;
}

void MetaDataVec::init(const std::vector<MDLabel> &labelsVector) {
    size_t col = 0;
    for (const auto label : labelsVector) {
        _label_to_col[label] = col;
        col++;
    }
    _no_columns = col;
}

int MetaDataVec::_labelIndex(MDLabel label) const {
    return this->_label_to_col[label];
}

const MDObject& MetaDataVec::_getObject(size_t i, MDLabel label) const {
    int labelIndex = _labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
    return this->_rows.at(i).at(labelIndex);
}

MDObject& MetaDataVec::_getObject(size_t i, MDLabel label) {
    int labelIndex = _labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
    return this->_rows.at(i).at(labelIndex);
}

int MetaDataVec::_rowIndex(size_t id) const {
    if (this->_id_to_index.find(id) == this->_id_to_index.end())
        return -1;
    return this->_id_to_index.at(id);
}

size_t MetaDataVec::_rowIndexSafe(size_t id) const {
    int i = this->_rowIndex(id);
    if (i == -1)
        throw ObjectDoesNotExist();
    return i;
}

void MetaDataVec::read(const FileName &inFile, const std::vector<MDLabel> *dediredLables, bool decomposeStack) {
    // FIXME: implement
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
        // FIXME: implement
        throw NotImplemented();
    } else {
        writeStar(_outFile, blockName, mode);
    }
}

void MetaDataVec::clear() {
    this->baseClear();
    this->_rows.clear();
    std::fill(this->_label_to_col.begin(), this->_label_to_col.end(), -1);
    this->_col_to_label.clear();
    this->_no_columns = 0;
}

void MetaDataVec::_setRow(const MDRow &row, size_t index) {
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
    // FIXME: should id be changed or kept same?
    MetaDataVecRow newRow;
    _rows.push_back(newRow);
    this->_setRow(row, _rows.size()-1);
    return getRowId(_rows.size()-1);
}

void MetaDataVec::addRows(const std::vector<MDRowVec> &rows) {
    for (const auto row : rows)
        this->addRow(row);
}

int MetaData::getMaxStringLength(const MDLabel thisLabel) const {
    return 0xFF;
}

bool MetaDataVec::setValueCol(const MDObject &mdValueIn) {
    int labelIndex = this->_labelIndex(mdValueIn.label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
    for (auto& vecRow : this->_rows)
        vecRow.at(labelIndex) = mdValueIn;
    return true;
}

bool MetaDataVec::setValue(const MDObject &mdValueIn, size_t id) {
    this->_getObject(this->_rowIndexSafe(id), mdValueIn.label) = mdValueIn;
    return true;
}

bool MetaDataVec::getValue(MDObject &mdValueOut, size_t id) const {
    mdValueOut = this->_getObject(this->_rowIndexSafe(id), mdValueOut.label);
    return true;
}

std::unique_ptr<MDRow> MetaDataVec::getRow(size_t id) {
    int i = this->_rowIndex(id);
    if (i < 0)
        return nullptr;
    return MemHelpers::make_unique<MDRowVec>(this->_rows[i], i, this->_label_to_col);
}

std::unique_ptr<const MDRow> MetaDataVec::getRow(size_t id) const {
    int i = this->_rowIndex(id);
    if (i < 0)
        return nullptr;
    return MemHelpers::make_unique<MDRowVec>(this->_rows[i], i, this->_label_to_col);
}

MDRowVec MetaDataVec::getRowVec(size_t id) {
    size_t i = this->_rowIndexSafe(id);
    return MDRowVec(this->_rows[i], i, this->_label_to_col);
}

const MDRowVec MetaDataVec::getRowVec(size_t id) const {
    size_t i = this->_rowIndexSafe(id);
    return MDRowVec(this->_rows[i], i, this->_label_to_col);
}


bool MetaDataVec::getRow(MDRow &row, size_t id) {
    size_t i = this->_rowIndexSafe(id);
    row = MDRowVec(this->_rows[i], i, this->_label_to_col);
    return true;
}

void MetaDataVec::getRow(MDRowVec &row, size_t id) {
    size_t i = this->_rowIndexSafe(id);
    row = MDRowVec(this->_rows[i], i, this->_label_to_col);
}

bool MetaDataVec::getRowValues(size_t id, std::vector<MDObject> &values) const {
    size_t i = this->_rowIndexSafe(id);
    values = this->_rows[i];
    return true;
}

size_t MetaDataVec::getRowId(size_t i) const {
    size_t id;
    this->_getObject(i, MDL_OBJID).getValue(id);
    return id;
}

size_t MetaDataVec::getRowId(const MetaDataVecRow& row) const {
    int labelIndex = _labelIndex(MDL_OBJID);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
    size_t id;
    row.at(labelIndex).getValue(id);
    return id;
}

void MetaDataVec::getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const {
    valuesOut.clear();
    int labelIndex = this->_labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
    for (const auto& vecRow : this->_rows)
        valuesOut.push_back(vecRow.at(labelIndex));
}

void MetaDataVec::setColumnValues(const std::vector<MDObject> &valuesIn) {
    for (size_t i = 0; i < std::min(valuesIn.size(), this->_rows.size()); i++) {
        int labelIndex = this->_labelIndex(valuesIn[i].label);
        if (labelIndex < 0)
            this->addLabel(valuesIn[i].label);
        labelIndex = this->_labelIndex(valuesIn[i].label);
        if (labelIndex < 0)
            throw ColumnDoesNotExist(); // internal error
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
        throw NotImplemented();
    if (this->_label_to_col[label] != -1)
        return true;

    this->_no_columns++;
    size_t column = this->_no_columns-1;
    this->_label_to_col[label] = column;
    this->_col_to_label[column] = label;
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
    return true;
}

bool MetaDataVec::keepLabels(const std::vector<MDLabel> &labels) {
    // FIXME: implement
    throw NotImplemented(); // not implemented yet
}

size_t MetaDataVec::addObject() {
    // FIXME: implement
    throw NotImplemented(); // not implemented yet
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

    this->_id_to_index.erase(i);
    this->_rows.erase(this->_rows.begin()+i);

    for (size_t j = i; j < this->_rows.size(); j++)
        this->_id_to_index[this->getRowId(i)]--;

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

void MetaDataVec::addItemId() {
    // FIXME: implement
    throw NotImplemented();
}

size_t MetaDataVec::firstRowId() const {
    return this->getRowId(0);
}

size_t MetaDataVec::lastRowId() const {
    return this->getRowId(this->size()-1);
}

bool MetaDataVec::_match(const MetaDataVecRow& row, const MDQuery& query) const {
    if (dynamic_cast<const MDValueRelational*>(&query) == nullptr)
        throw NotImplemented(); // MDValueRange, MDExpression, MDMultiQuery not implemented yet
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
            objectsOut.push_back(this->getRowId(row));
}

void MetaDataVec::findObjects(std::vector<size_t> &objectsOut, int limit) const {
    objectsOut.clear();
    for (size_t i = 0; i < std::min<size_t>(limit, this->size()); i++)
        objectsOut.push_back(this->getRowId(this->_rows[i]));
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

/*
void _writeRows(std::ostream &os) const;*/

void MetaDataVec::writeStar(const FileName &outFile, const String & blockName, WriteModeMetaData mode) const {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::write(std::ostream &os, const String &blockName, WriteModeMetaData mode) const {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::print() const {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::append(const FileName &outFile) const {
    // TODO
    throw NotImplemented();
}

/*bool existsBlock(const FileName &_inFile);
void readStar(const FileName &inFile,
              const std::vector<MDLabel> *desiredLabels = nullptr,
              const String & blockName=DEFAULT_BLOCK_NAME,
              bool decomposeStack=true);
void readXML(const FileName &inFile,
             const std::vector<MDLabel> *desiredLabels= nullptr,
             const String & blockRegExp=DEFAULT_BLOCK_NAME,
             bool decomposeStack=true);

void read(const FileName &inFile, const std::vector<MDLabel> *desiredLabels = nullptr, bool decomposeStack = true) override;
void readPlain(const FileName &inFile, const String &labelsString, const String &separator = " ");
void addPlain(const FileName &inFile, const String &labelsString, const String &separator=" ");
*/

double MetaDataVec::getColumnMax(MDLabel column) {
    // TODO
    throw NotImplemented();
}

double MetaDataVec::getColumnMin(MDLabel column) {
    // TODO
    throw NotImplemented();
}

/*
void operate(const String &expression);
void replace(const MDLabel label, const String &oldStr, const String &newStr);*/

void MetaDataVec::randomize(const MetaData &MDin) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::removeDuplicates(MetaData &MDin, MDLabel label) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::removeDisabled() {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::sort(MetaDataVec &MDin, const MDLabel sortLabel, bool asc, int limit, int offset) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::sort(MetaDataVec &MDin, const String &sortLabel, bool asc, int limit, int offset) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::split(size_t n, std::vector<MetaDataVec> &results, const MDLabel sortLabel) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::selectSplitPart(const MetaData &mdIn, size_t n, size_t part,
                                   const MDLabel sortLabel) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                             const MDLabel sortLabel) {
    // TODO
    throw NotImplemented();
}

/*void makeAbsPath(const MDLabel label=MDL_IMAGE);*/


void MetaDataVec::fillExpand(MDLabel label) {
    // FIXME
    throw NotImplemented();
}

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
        row.push_back(MDObject(this->_col_to_label.at(i)));
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

        if (static_cast<size_t>(labelsrci) < row.size())
            row[labeldesti] = row[labelsrci];
        // else row[labeldesti] is empty MDObject (from previous if)
    }
}

void MetaDataVec::copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::renameColumn(MDLabel oldLabel, MDLabel newLabel) {
    // TODO
    throw NotImplemented();
}

void MetaDataVec::renameColumn(const std::vector<MDLabel> &oldLabel,
                               const std::vector<MDLabel> &newLabel) {
    // TODO
    throw NotImplemented();
}
