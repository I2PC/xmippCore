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

void MetaDataVec::read(const FileName &inFile, const std::vector<MDLabel> *dediredLables, bool decomposeStack) {
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

void MetaDataVec::clear() {
    this->baseClear();
    this->_rows.clear();
    std::fill(this->_label_to_col.begin(), this->_label_to_col.end(), -1);
    this->_col_to_label.clear();
    this->_no_columns = 0;
}

template <typename T>
size_t MetaDataVec::addRow(const T &row) {
    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(labeli);
        if (row.containsLabel(label) && !this->containsLabel(label))
            this->addLabel(label);
    }

    MetaDataVecRow newRow;

    for (size_t coli = 0; coli < this->_no_columns; coli++)
        newRow.push_back({this->_col_to_label[coli]});

    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; ++labeli) {
        MDLabel label = static_cast<MDLabel>(labeli);
        if (row.containsLabel(label)) {
            size_t ourCol = _label_to_col[label];
            row.getValue(newRow[ourCol]);
        }
    }

    _rows.push_back(newRow);
    return getRowId(_rows.size()-1); // FIXME: should id be changed or kept same?
}

void MetaDataVec::addRows(const std::vector<MDRowVec> &rows) {
    for (const auto row : rows)
        this->addRow(row);
}

int MetaData::getMaxStringLength(const MDLabel thisLabel) const {
    return 0xFF;
}


/*bool MetaData::setValueCol(const MDObject &mdValueIn) {
}*/


/*bool setValue(const MDObject &mdValueIn, size_t id);

bool getValue(MDObject &mdValueOut, size_t id) const;

std::unique_ptr<MDRow> getRow(size_t id) const override;
bool getRow(MDRow &row, size_t id) const override;
bool getRow(MDRowVec &row, size_t id) const;
bool getRowValues(size_t id, std::vector<MDObject> &values) const override;
size_t getRowId(size_t i) const;
void getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const override;
void setColumnValues(const std::vector<MDObject> &valuesIn) override;*/

bool MetaDataVec::setRow(const MDRow &row, size_t id) {
}

bool MetaDataVec::setRow(const MDRowConst &row, size_t id) {
}

bool MetaDataVec::setRow(const MDRowVec &row, size_t id) {
}

bool MetaDataVec::setRow(const MDRowVecConst &row, size_t id) {
}

/*bool setValueFromStr(const MDLabel label, const String &value, size_t id);
bool getStrFromValue(const MDLabel label, String &strOut, size_t id) const;
bool isEmpty() const override;
size_t size() const override;
bool containsLabel(const MDLabel label) const override;
bool addLabel(const MDLabel label, int pos = -1) override;
bool removeLabel(const MDLabel label) override;
bool keepLabels(const std::vector<MDLabel> &labels) override;

size_t addObject() override;

void importObject(const MetaData &md, const size_t id, bool doClear=true) override;
void importObjects(const MetaData &md, const std::vector<size_t> &objectsToAdd, bool doClear=true) override;
void importObjects(const MetaData &md, const MDQuery &query, bool doClear=true) override;

bool removeObject(size_t id) override;

void removeObjects(const std::vector<size_t> &toRemove) override;

int removeObjects() override;
int removeObjects(const MDQuery&) override;

void addItemId();

size_t firstRowId() const override;
size_t firstObject(const MDQuery&) const override;
size_t lastRowId() const override;

void findObjects(std::vector<size_t> &objectsOut, const MDQuery &query) const override;
void findObjects(std::vector<size_t> &objectsOut, int limit = -1) const override;

size_t countObjects(const MDQuery&) const override;
bool containsObject(size_t objectId) const override;
bool containsObject(const MDQuery&) const override;

bool containsObject(size_t objectId);

void _writeRows(std::ostream &os) const;

void writeStar(const FileName &outFile,const String & blockName="", WriteModeMetaData mode=MD_OVERWRITE) const;
void write(const FileName &outFile, WriteModeMetaData mode=MD_OVERWRITE) const;
void write(std::ostream &os, const String & blockName="",WriteModeMetaData mode=MD_OVERWRITE) const;
void print() const;
void append(const FileName &outFile) const;
bool existsBlock(const FileName &_inFile);
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


double getColumnMax(MDLabel column);
double getColumnMin(MDLabel column);


void operate(const String &expression);
void replace(const MDLabel label, const String &oldStr, const String &newStr);
void randomize(const MetaData &MDin);
void removeDuplicates(MetaData &MDin, MDLabel label=MDL_UNDEFINED);
void removeDisabled();

void sort(MetaDataVec &MDin,
          const MDLabel sortLabel,
          bool asc=true,
          int limit=-1,
          int offset=0);

void sort(MetaDataVec &MDin, const String &sortLabel, bool asc=true, int limit=-1, int offset=0);

void split(size_t n, std::vector<MetaDataVec> &results,
           const MDLabel sortLabel=MDL_OBJID);

void selectSplitPart(const MetaData &mdIn,
                     size_t n, size_t part,
                     const MDLabel sortLabel=MDL_OBJID) override;

void selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel=MDL_OBJID) override;

void selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                const MDLabel sortLabel=MDL_OBJID) override;

void makeAbsPath(const MDLabel label=MDL_IMAGE);


void fillExpand(MDLabel label) override;
void fillConstant(MDLabel label, const String &value) override;
void fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3=0.) override;
void fillLinear(MDLabel label, double initial, double step) override;

void copyColumn(MDLabel labelDest, MDLabel labelSrc) override;
void copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc) override;

void renameColumn(MDLabel oldLabel, MDLabel newLabel) override;
void renameColumn(const std::vector<MDLabel> &oldLabel,
        const std::vector<MDLabel> &newLabel) override;*/
