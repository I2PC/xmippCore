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
#include <cassert>
#include <fstream>
#include "metadata_vec.h"
#include "metadata_generator.h"
#include "xmipp_image.h"

#include <sys/stat.h>
#include <fcntl.h>
#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif

#define END_OF_LINE() ((char*) memchr (iter, '\n', end-iter))

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
    this->clear();
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
    return this->_getObject(this->_rows.at(i), label);
}

MDObject& MetaDataVec::_getObject(size_t i, MDLabel label) {
    return this->_getObject(this->_rows.at(i), label);
}

const MDObject& MetaDataVec::_getObject(const MetaDataVecRow& row, MDLabel label) const {
    int labelIndex = this->_labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
    return row.at(labelIndex);
}

MDObject& MetaDataVec::_getObject(MetaDataVecRow& row, MDLabel label) const {
    int labelIndex = this->_labelIndex(label);
    if (labelIndex < 0)
        throw ColumnDoesNotExist();
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
        throw ObjectDoesNotExist();
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
    this->_isColumnFormat = true;

    if (extFile == "xml")
        this->readXML(_inFile, desiredLabels, blockName, decomposeStack);
    else if (extFile == "sqlite")
        throw NotImplemented("Reading from .sqlite file into MetaDataVec not implemented!");
    else
        this->readStar(filename, desiredLabels, blockName, decomposeStack);

    this->eFilename = filename;
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
    this->baseClear();
    this->_rows.clear();
    std::fill(this->_label_to_col.begin(), this->_label_to_col.end(), -1);
    this->_col_to_label.clear();
    this->_id_to_index.clear();
    this->_no_columns = 0;
    this->_next_id = 0;
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

    if (!row.containsLabel(MDL_OBJID)) {
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

    if (rowId == this->_next_id)
        this->_next_id++;
    this->_id_to_index[rowId] = _rows.size()-1;
    return rowId;
}

void MetaDataVec::addRows(const std::vector<MDRowVec> &rows) {
    for (const auto row : rows)
        this->addRow(row);
}

int MetaDataVec::getMaxStringLength(const MDLabel thisLabel) const {
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
    if (!this->containsLabel(mdValueIn.label))
        this->addLabel(mdValueIn.label);

    MetaDataVecRow& row = this->_rows[this->_rowIndexSafe(id)];
    this->_expand(row, mdValueIn.label);

    row[this->_labelIndex(mdValueIn.label)] = mdValueIn;
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
        throw NotImplemented("addLabel to -1 not implemented");
    if (this->_label_to_col[label] != -1)
        return true;

    this->_no_columns++;
    size_t column = this->_no_columns-1;
    this->_label_to_col[label] = column;
    this->_col_to_label.push_back(label);
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
    throw NotImplemented("keepLabels not implemented"); // not implemented yet
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
    throw NotImplemented("addItemId not implemented");
}

size_t MetaDataVec::firstRowId() const {
    return this->getRowId(0);
}

size_t MetaDataVec::lastRowId() const {
    return this->getRowId(this->size()-1);
}

bool MetaDataVec::_match(const MetaDataVecRow& row, const MDQuery& query) const {
    if (dynamic_cast<const MDValueRelational*>(&query) == nullptr)
        throw NotImplemented("_match not implemented"); // MDValueRange, MDExpression, MDMultiQuery not implemented yet
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


void MetaDataVec::_writeRows(std::ostream &os) const {
    for (const MetaDataVecRow& row : this->_rows) {
        for (size_t i = 0; i < MDL_LAST_LABEL; i++) {
            const MDLabel label = static_cast<MDLabel>(i);
            if ((label != MDL_STAR_COMMENT) && (this->_label_to_col[i] > -1) &&
                (this->_labelIndex(label) < static_cast<int>(this->_rows[0].size()))) {
                os.width(1);
                this->_getObject(row, label).toStream(os);
                os << " ";
            }
        }
        os << '\n';
    }
}

void MetaDataVec::writeStar(const FileName &outFile, const String & blockName, WriteModeMetaData mode) const {
    // FIXME: this method is compeletely same as in MetaDataDb
    // Move to MetaData?
#ifdef XMIPP_MMAP
    if (outFile.hasImageExtension())
        REPORT_ERROR(ERR_IO,"MetaData:writeStar Trying to write metadata with image extension");

    struct stat file_status;
    int fd;
    char *map = NULL;
    char *tailMetadataFile = NULL; // auxiliary variable to keep metadata file tail in memory
    size_t size=-1;
    char *target, * target2 = NULL;

    //check if file exists or not block name has been given
    //in our format no two identical data_xxx strings may exists
    if (mode == MD_APPEND) {
        if (blockName.empty() || !outFile.exists())
            mode = MD_OVERWRITE;
        else {
            //does blockname exists?
            //remove it from file in this case
            // get length of file:
            if(stat(outFile.c_str(), &file_status) != 0)
                REPORT_ERROR(ERR_IO_NOPATH,"MetaData:writeStar can not get filesize for file "+outFile);
            size = file_status.st_size;
            if (size == 0)
                mode = MD_OVERWRITE;
        }

        if (mode == MD_APPEND) { //size=0 for /dev/stderr
            fd = open(outFile.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
            if (fd == -1)
                REPORT_ERROR(ERR_IO_NOPATH,"MetaData:writeStar can not read file named "+outFile);

            map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (map == MAP_FAILED)
                REPORT_ERROR(ERR_MEM_BADREQUEST,"MetaData:writeStar can not map memory ");

            // Is this a metadata formatted FILE
            if(strncmp(map,FileNameVersion.c_str(),FileNameVersion.length()) != 0) {
                mode = MD_OVERWRITE;
            } else {
                //block name
                String _szBlockName = formatString("\ndata_%s\n", blockName.c_str());
                size_t blockNameSize = _szBlockName.size();

                //search for the string
                target = (char *) _memmem(map, size, _szBlockName.c_str(), blockNameSize);

                if (target != NULL)
                {
                    target2 = (char *) _memmem(target+1, size - (target - map), "\ndata_", 6);

                    if (target2 != NULL)
                    {
                        //target block is not the last one, so we need to
                        //copy file from target2 to auxiliary memory and truncate
                        tailMetadataFile = (char *) malloc( ((map + size) - target2));
                        memmove(tailMetadataFile,target2, (map + size) - target2);
                    }
                    if (ftruncate(fd, target - map+1)==-1)  //truncate rest of the file
                        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot truncate file");
                }
            }
            close(fd);

            if (munmap(map, size) == -1)
                REPORT_ERROR(ERR_MEM_NOTDEALLOC, "MetaData:writeStar, Can not unmap memory");
        }
    }

    std::ios_base::openmode openMode = (mode == MD_OVERWRITE) ? std::ios_base::out : std::ios_base::app;
    std::ofstream ofs(outFile.c_str(), openMode);

    write(ofs, blockName, mode);

    if (tailMetadataFile != NULL) {
        //append memory buffer to file
        //may a cat a buffer to a ofstream
        ofs.write(tailMetadataFile,(map + size) - target2);
        free(tailMetadataFile);
    }
    ofs.close();

#else

    REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif
}

void MetaDataVec::write(std::ostream &os, const String &blockName, WriteModeMetaData mode) const {
    if (mode==MD_OVERWRITE)
        os << FileNameVersion << " * "// << (isColumnFormat ? "column" : "row")
        << '\n' // write which type of format (column or row) and the path;
        << WordWrap(this->_comment, line_max);     // write md comment in the 2nd comment line of header

    // write data block
    String _szBlockName("data_");
    _szBlockName += blockName;

    if (this->_isColumnFormat) {
        // write md columns in 3rd comment line of the header
        os << _szBlockName << '\n';
        os << "loop_" << '\n';
        for (size_t i = 0; i < MDL_LAST_LABEL; i++)
            if ((i != MDL_STAR_COMMENT) && (this->_label_to_col[i] > -1))
                os << " _" << MDL::label2Str(static_cast<MDLabel>(i)) << '\n';
        _writeRows(os);

        //Put the activeObject to the first, if exists
    } else { // row format
        os << _szBlockName << '\n';

        // Print single object
        assert(this->_rows.size() == 1);

        for (size_t i = 0; i < MDL_LAST_LABEL; i++) {
            const MDLabel label = static_cast<MDLabel>(i);
            if ((label != MDL_STAR_COMMENT) && (this->_label_to_col[i] > -1) &&
                (this->_labelIndex(label) < static_cast<int>(this->_rows[0].size()))) {
                this->_getObject(0, label).toStream(os);
                os << '\n';
            }
        }
    }
}

void MetaDataVec::append(const FileName &outFile) const {
    // TODO
    throw NotImplemented("append not implemented");
}

bool MetaDataVec::existsBlock(const FileName &_inFile) {
    // TODO
    throw NotImplemented("existsBlock not implemented");
}

void MetaDataVec::readStar(const FileName &filename, const std::vector<MDLabel> *desiredLabels,
                           const String &blockRegExp, bool decomposeStack) {
    // FIXME: move this to MetaData?
    // First try to open the file as a metadata
    size_t id;
    FileName inFile = filename.removeBlockName();

    if (!(isMetadataFile = inFile.isMetaData())) { //if not a metadata, try to read as image or stack
        Image<char> image;
        if (decomposeStack) // If not decomposeStack it is no necessary to read the image header
            image.read(filename, HEADER);
        if (!decomposeStack || image().ndim == 1) { // single image; !decomposeStack must be first
            id = this->addObject();
            MetaData::setValue(MDL_IMAGE, filename, id);
            MetaData::setValue(MDL_ENABLED, 1, id);
        } else { // stack
            FileName fnTemp;
            for (size_t i = 1; i <= image().ndim; ++i) {
                fnTemp.compose(i, filename);
                id = addObject();
                MetaData::setValue(MDL_IMAGE, fnTemp, id);
                MetaData::setValue(MDL_ENABLED, 1, id);
            }
        }
        return;
    }

    std::ifstream is(inFile.c_str(), std::ios_base::in);
    std::stringstream ss;
    String line, token,_comment;
    std::vector<MDObject*> columnValues;

    getline(is, line); //get first line to identify the type of file

    if (is.fail())
        REPORT_ERROR(ERR_IO_NOTEXIST, formatString("MetaDataDb::read: File doesn't exists: %s", inFile.c_str()) );

    bool useCommentAsImage = false;
    this->_inFile = inFile;
    bool oldFormat = true;

    is.seekg(0, std::ios::beg);//reset the stream position to the beginning to start parsing

    if (line.find(FileNameVersion) != String::npos ||
        eFilename.getExtension() == "xmd" || eFilename.getExtension() == "star") {
        oldFormat = false;
        _comment.clear();

        // Skip comment parsing if we found the data key in the first line
        if (line.find("data_") != 0) {
            // Read comment
            //        is.ignore(256,'#');//format line
            is.ignore(256, '\n');//skip first line
            bool addspace = false;
            while (1) {
                getline(is, line);
                trim(line);
                if (line[0] == '#')
                {
                    line[0] = ' ';
                    trim(line);
                    if (addspace)
                        _comment += " " + line;
                    else
                        _comment += line;
                    addspace = true;
                } else
                    break;
            }
            this->setComment(_comment);
        }

        //map file
        int fd;
        BUFFER_CREATE(bufferMap);
        mapFile(inFile, bufferMap.begin, bufferMap.size, fd);

        BLOCK_CREATE(block);
        regex_t re;
        int rc = regcomp(&re, (blockRegExp+"$").c_str(), REG_EXTENDED|REG_NOSUB);
        if (blockRegExp.size() && rc != 0)
            REPORT_ERROR(ERR_ARG_INCORRECT, formatString("Pattern '%s' cannot be parsed: %s",
                         blockRegExp.c_str(), inFile.c_str()));
        BUFFER_COPY(bufferMap, buffer);
        bool firstBlock = true;
        bool singleBlock = blockRegExp.find_first_of(".[*+")==String::npos;

        String blockName;

        while (nextBlock(buffer, block)) {
            //startingPoint, remainingSize, firstData, secondData, firstloop))
            BLOCK_NAME(block, blockName);
            if (blockRegExp.size() == 0 || regexec(&re, blockName.c_str(), (size_t) 0, NULL, 0) == 0) {
                //Read column labels from the datablock that starts at firstData
                //Label ends at firstloop
                if ((_isColumnFormat = (block.loop != NULL))) {
                    _readColumnsStar(block, columnValues, desiredLabels, firstBlock);
                    // If block is empty, makes block.loop and block.end equal
                    if (block.loop == (block.end + 1))
                        block.loop--;
                    _readRowsStar(block, columnValues, desiredLabels);
                } else {
                    id = this->addObject();
                    _parsedLines = 1;
                    this->_readColumnsStar(block, columnValues, desiredLabels, firstBlock, id);
                }
                firstBlock = false;

                if (singleBlock)
                    break;
            }
        }

        unmapFile(bufferMap.begin, bufferMap.size, fd);
        regfree(&re);
        if (firstBlock)
            REPORT_ERROR(ERR_MD_BADBLOCK, formatString("Block: '%s': %s",
                         blockRegExp.c_str(), inFile.c_str()));

    } else if (line.find("Headerinfo columns:") != String::npos) {
        // This looks like an old DocFile, parse header
        std::cerr << "WARNING: ** You are using an old file format (DOCFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        is.ignore(256, ':'); //ignore all until ':' to start parsing column labels
        getline(is, line);
        ss.str(line);
        columnValues.push_back(new MDObject(MDL_UNDEFINED));
        columnValues.push_back(new MDObject(MDL_UNDEFINED));

        this->addLabel(MDL_IMAGE);
        this->_readColumns(ss, columnValues, desiredLabels);
        useCommentAsImage = true;

    } else {
        std::cerr << "WARNING: ** You are using an old file format (SELFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        // I will assume that is an old SelFile, so only need to add two columns
        columnValues.push_back(new MDObject(MDL_IMAGE));//addLabel(MDL_IMAGE);
        columnValues.push_back(new MDObject(MDL_ENABLED));//addLabel(MDL_ENABLED);
    }

    if (oldFormat)
        this->_readRows(is, columnValues, useCommentAsImage);

    //free memory of column values
    int nCols = columnValues.size();
    for (int i = 0; i < nCols; ++i)
        delete columnValues[i];

    is.close();
}

void MetaDataVec::_readRows(std::istream& is, std::vector<MDObject*>& columnValues, bool useCommentAsImage) {
    // FIXME: move this to MetaData?
    String line = "";
    while (!is.eof() && !is.fail()) {
        // Move until the ';' or the first alphanumeric character
        while (is.peek() != ';' && isspace(is.peek()) && !is.eof())
            is.ignore(1);
        if (!is.eof()) {
            if (is.peek() == ';') { //is a comment
                is.ignore(1); //ignore the ';'
                getline(is, line);
                trim(line);
            } else if (!isspace(is.peek())) {
                size_t id = addObject();
                if (line != "") { //this is for old format files
                    if (!useCommentAsImage)
                        this->setValue(MDObject(MDL_STAR_COMMENT, line), id);
                    else
                        this->setValue(MDObject(MDL_IMAGE, line), id);
                }
                int nCol = columnValues.size();
                for (int i = 0; i < nCol; ++i)
                    this->_parseObject(is, *(columnValues[i]), id);
            }
        }
    }
}

/* This function parses rows data in START format
 */
void MetaDataVec::_readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues,
                                const std::vector<MDLabel> *desiredLabels) {
    // FIXME: move this to MetaData?
    String line;
    std::stringstream ss;
    size_t n = block.end - block.loop;
    bool firstTime = true;

    if (n == 0)
        return;

    char * buffer = new char[n];
    memcpy(buffer, block.loop, n);
    char *iter = buffer, *end = iter + n, * newline = NULL;
    _parsedLines = 0; //Check how many lines the md have

    while (iter < end) { //while there are data lines
        //Adding \n position and check if NULL at the same time
        if (!(newline = END_OF_LINE()))
            newline = end;
        line.assign(iter, newline - iter);
        trim(line);

        if (!line.empty() && line[0] != '#') {
            //_maxRows would be > 0 if we only want to read some
            // rows from the md for performance reasons...
            // anyway the number of lines will be counted in _parsedLines
            if (_maxRows == 0 || _parsedLines < _maxRows) {
                std::stringstream ss(line);
                this->_parseObjects(ss, columnValues, desiredLabels, firstTime);
                firstTime = false;
            }
            _parsedLines++;
        }
        iter = newline + 1; //go to next line
    }

    delete[] buffer;
}

void MetaDataVec::_parseObjects(std::istream &is, std::vector<MDObject*> &columnValues,
                                const std::vector<MDLabel> *desiredLabels, bool firstTime) {
    // FIXME: move this to MetaData?
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

/* Helper function to parse an MDObject and set its value.
 * The parsing will be from an input stream(istream)
 * and if parsing fails, an error will be raised
 */
void MetaDataVec::_parseObject(std::istream &is, MDObject &object, size_t id) {
    // FIXME: move this to MetaData?
    object.fromStream(is);
    if (is.fail()) {
       String errorMsg = formatString("MetaData: Error parsing column '%s' value.", MDL::label2Str(object.label).c_str());
       object.failed = true;
       std::cerr << "WARNING: " << errorMsg << std::endl;
       //REPORT_ERROR(ERR_MD_BADLABEL, (String)"read: Error parsing data column, expecting " + MDL::label2Str(object.label));
    } else {
        if (object.label != MDL_UNDEFINED)
            this->setValue(object, id);
    }
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
    // TODO
    throw NotImplemented("randomize not implemented");
}

void MetaDataVec::removeDuplicates(MetaData &MDin, MDLabel label) {
    // TODO
    throw NotImplemented("removeDuplicates not implemented");
}

void MetaDataVec::removeDisabled() {
    // TODO
    throw NotImplemented("removeDisabled not implemented");
}

void MetaDataVec::sort(MetaDataVec &MDin, const MDLabel sortLabel, bool asc, int limit, int offset) {
    // TODO
    throw NotImplemented("sort not implemented");
}

void MetaDataVec::sort(MetaDataVec &MDin, const String &sortLabel, bool asc, int limit, int offset) {
    // TODO
    throw NotImplemented("sort not implemented");
}

void MetaDataVec::split(size_t n, std::vector<MetaDataVec> &results, const MDLabel sortLabel) {
    // TODO
    throw NotImplemented("split not implemented");
}

void MetaDataVec::selectSplitPart(const MetaData &mdIn, size_t n, size_t part,
                                   const MDLabel sortLabel) {
    // TODO
    throw NotImplemented("selectSplitPart not implemented");
}

void MetaDataVec::selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel) {
    // TODO
    throw NotImplemented("selectRandomSubset not implemented");
}

void MetaDataVec::selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                             const MDLabel sortLabel) {
    // TODO
    throw NotImplemented("selectPart not implemented");
}

/*void makeAbsPath(const MDLabel label=MDL_IMAGE);*/


void MetaDataVec::fillExpand(MDLabel label) {
    // FIXME
    throw NotImplemented("fillExpand not implemented");
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
    throw NotImplemented("copyColumnTo not implemented");
}

void MetaDataVec::renameColumn(MDLabel oldLabel, MDLabel newLabel) {
    // TODO
    throw NotImplemented("removeColumn not implemented");
}

void MetaDataVec::renameColumn(const std::vector<MDLabel> &oldLabel,
                               const std::vector<MDLabel> &newLabel) {
    // TODO
    throw NotImplemented("renameColumn not implemented");
}


bool MetaDataVec::operator==(const MetaDataVec& op) const {
    // This comparison ignores order of labels, everything else must be same,
    // even order of rows and ids.

    if (this->_rows.size() != op._rows.size())
        return false;

    for (size_t labeli = 0; labeli < MDL_LAST_LABEL; labeli++)
        if ((this->_label_to_col[labeli] > -1) != (op._label_to_col[labeli] > -1))
            return false;

    for (size_t i = 0; i < this->_rows.size(); i++) {
        for (size_t labeli = 0; labeli < MDL_LAST_LABEL; labeli++) {
            int thisLabelRowI = this->_label_to_col[labeli];
            int opLabelRowI = op._label_to_col[labeli];
            if (thisLabelRowI > -1) {
                if ((static_cast<size_t>(thisLabelRowI) < this->_rows[i].size()) !=
                    (static_cast<size_t>(opLabelRowI) < op._rows[i].size()))
                    return false; // item present in one row, but not other
                if (static_cast<size_t>(thisLabelRowI) >= this->_rows[i].size())
                    continue; // label not present in both rows
                if (!this->_rows[i][thisLabelRowI].eq(op._rows[i][opLabelRowI], this->precision()))
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
            out.push_back(static_cast<MDLabel>(i));
    return out;
}

std::ostream& operator<<(std::ostream& o, const MetaDataVec& md) {
    md.write(o);
    return o;
}
