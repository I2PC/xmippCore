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

#include <fstream>
#include <random>
#include <algorithm>
#include <cassert>

#include "metadata_db.h"
#include "xmipp_image.h"
#include "metadata_sql.h"
#include "metadata_generator.h"
#include "xmipp_funcs.h"


//-----Constructors and related functions ------------
void MetaDataDb::_clear(bool onlyData)
{
    if (onlyData)
    {
        myMDSql->deleteObjects();
    }
    else
    {
        MetaData::clear();
        _activeLabels.clear();
        myMDSql->clearMd();
    }
}//close clear

void MetaDataDb::clear()
{
    init({});
}

void MetaDataDb::init(const std::vector<MDLabel> &labelsVector)
{
    _clear();
    _maxRows = 0; //by default read all rows
    _parsedLines = 0; //no parsed line;
    _activeLabels = labelsVector;
    //Create table in database
    myMDSql->createMd();
    _precision = 100;
    isMetadataFile = false;
}//close init

void MetaDataDb::copyMetadata(const MetaDataDb &md, bool copyObjects)
{
    if (this == &md) //not sense to copy same metadata
        return;
    init(md._activeLabels);
    copyInfo(md);
    if (!md._activeLabels.empty())
    {
        if (copyObjects)
            md.myMDSql->copyObjects(this);
    }
    else
    {
        int n = md.size();
        for (int i = 0; i < n; i++)
            addObject();
    }
}

bool MetaDataDb::setValue(const MDObject &mdValueIn, size_t id)
{
    if (id == BAD_OBJID)
    {
        REPORT_ERROR(ERR_MD_NOACTIVE, "setValue: please provide objId other than -1");
        exit(1);
    }
    //add label if not exists, this is checked in addlabel
    addLabel(mdValueIn.label);
    return myMDSql->setObjectValue(id, mdValueIn);
}

bool MetaDataDb::setValueCol(const MDObject &mdValueIn)
{
    //add label if not exists, this is checked in addlabel
    addLabel(mdValueIn.label);
    return myMDSql->setObjectValue(mdValueIn);
}

bool MetaDataDb::getValue(MDObject &mdValueOut, size_t id) const
{
    if (!containsLabel(mdValueOut.label))
        return false;

    if (id == BAD_OBJID)
        REPORT_ERROR(ERR_MD_NOACTIVE, "getValue: please provide objId other than -1");

    return myMDSql->getObjectValue(id, mdValueOut);
}

void MetaDataDb::getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const
{
    MDObject mdValueOut(label);
    std::vector<size_t> objectsId;
    findObjects(objectsId);
    size_t n = objectsId.size();
    valuesOut.resize(n,mdValueOut);
    for (size_t i = 0; i < n; ++i)
    {
        getValue(mdValueOut, objectsId[i]);
        valuesOut[i] = mdValueOut;
    }
}

void MetaDataDb::setColumnValues(const std::vector<MDObject> &valuesIn)
{
    bool addObjects=false;
    if (size()==0)
        addObjects=true;
    if (valuesIn.size()!=size() && !addObjects)
        REPORT_ERROR(ERR_MD_OBJECTNUMBER,"Input vector must be of the same size as the metadata");
    if (!addObjects)
    {
        size_t n = 0;
        for (size_t objId : this->ids())
            setValue(valuesIn[n++], objId);
    }
    else
    {
        size_t nmax=valuesIn.size();
        for (size_t n=0; n<nmax; ++n)
            setValue(valuesIn[n],addObject());
    }
}

bool MetaDataDb::bindValue(size_t id) const
{
    bool success=true;

    // Prepare statement.
    if (!myMDSql->bindStatement( id))
    {
        success = false;
    }

    return success;
}

bool MetaDataDb::initGetRow(bool addWhereClause) const
{
    bool success=true;

    // Prepare statement.
    if (!myMDSql->initializeSelect( addWhereClause, this->_activeLabels))
    {
        success = false;
    }

    return success;
}

bool MetaDataDb::execGetRow(MDRow &row) const
{
    std::vector<MDObject> mdValues;
    mdValues.reserve(this->_activeLabels.size());

    row.clear();

    bool success = myMDSql->getObjectsValues(this->_activeLabels, mdValues);
    if (success)
        for (const auto &obj : mdValues)
            row.setValue(obj);

    return success;
}

void MetaDataDb::finalizeGetRow(void) const
{
    myMDSql->finalizePreparedStmt();
}

std::vector<MDObject> MetaDataDb::getObjectsForActiveLabels() const {
    // get active labels
    std::vector<MDObject> values;
    const auto &labels = this->_activeLabels;
    values.reserve(labels.size());
    for (auto &l : labels) {
        values.emplace_back(l);
    }
    return values;
}

bool MetaDataDb::getAllRows(std::vector<MDRowSql> &rows) const
{
    std::vector<std::vector<MDObject>> rawRows;
    rawRows.reserve(this->size());
    auto columns = getObjectsForActiveLabels();
    if ( !  sqlUtils::select(myMDSql->db,
            myMDSql->tableName(myMDSql->tableId),
            columns,
            rawRows)) return false;

    rows.clear();
    const auto noOfRows = rawRows.size();
    rows.resize(noOfRows);
    for (size_t i = 0; i < noOfRows; ++i) {
        auto &row = rows.at(i);
        const auto &vals = rawRows.at(i);
        // fill the row
        for (auto &v : vals) {
            row.setValue(v);
        }
    }
    return true;
}

std::unique_ptr<MDRow> MetaDataDb::getRow(size_t id) {
    std::unique_ptr<MDRowSql> row(new MDRowSql());
    if (!getRow(*row, id))
        return nullptr;
    return std::move(row);
}

std::unique_ptr<const MDRow> MetaDataDb::getRow(size_t id) const {
    std::unique_ptr<MDRowSql> row(new MDRowSql());
    if (!getRow(*row, id))
        return nullptr;
    return std::move(row);
}

MDRowSql MetaDataDb::getRowSql(size_t id) {
    MDRowSql row;
    if (!getRow(row, id))
        throw ObjectDoesNotExist(id, getFilename());
    return row;
}

const MDRowSql MetaDataDb::getRowSql(size_t id) const {
    MDRowSql row;
    if (!getRow(row, id))
        throw ObjectDoesNotExist(id, getFilename());
    return row;
}

bool MetaDataDb::getRow(MDRowSql &row, size_t id) const
{
    if (id == BAD_OBJID)
        REPORT_ERROR(ERR_MD_NOACTIVE, "getValue: please provide objId other than -1");
    // clear whatever is there now
    row.clear();
    // get active labels
    auto values = getObjectsForActiveLabels();
    // get values from the row
    if ( ! sqlUtils::select(id,
            myMDSql->db,
            myMDSql->tableName(myMDSql->tableId),
            values)) return false;
    // fill them
    for (auto &v : values)
        row.setValue(v);
    return true;
}

bool MetaDataDb::getRow2(MDRow &row, size_t id) const
{
    bool success=true;

    // Clear row.
    row.clear();

    // Initialize SELECT.
    success = this->initGetRow( true);
    if (success)
    {
        bindValue( id);

        // Execute SELECT.
        success = execGetRow( row);

        // Finalize SELECT.
        finalizeGetRow();
    }

    return(success);
}

bool MetaDataDb::setRow(const MDRow &row, size_t id) 
{
    if (row.empty()) {
        return true;
    }
    addMissingLabels(row);

    // create mask of valid labels
    std::vector<MDLabel> labels;
    labels.reserve(row.size());
    for (const MDObject* obj : row) {
        labels.emplace_back(obj->label);
    }
    // extract values to be added
    std::vector<const MDObject*> vals;
    vals.reserve(row.size());
    for (const auto &l : labels) {
        vals.emplace_back(row.getObject(l));
    }
    // update values to db
    return sqlUtils::update(vals, MDSql::db,
                    myMDSql->tableName(myMDSql->tableId), id);
}


bool MetaDataDb::initAddRow(const MDRow &row)
{
    int     j=0;                    // Loop counter.
    bool    success=true;               // Return value.
    std::vector<MDLabel>    labels;     // Columns labels.
    std::vector<MDObject*>  mdValues;   // Vector to store values.

    // Set vector size.
    labels.resize(row.size());

    // Get labels.
    j=0;
    for (const MDObject* obj : row) {
        addLabel(obj->label);
        labels[j] = obj->label;
        j++;
    }
    labels.resize(j);

    // Prepare statement (mdValues is not used).
    if (!myMDSql->initializeInsert( &labels, mdValues))
    {
        std::cerr << "initAddRow: error executing myMDSql->initializeInsert" << std::endl;
        success = false;
    }

    return success;
}


bool MetaDataDb::execAddRow(const MDRow &row)
{
    int j = 0;
    bool success = true;
    std::vector<const MDObject*> mdValues;

    // Set values vector size.
    mdValues.resize(row.size());

    // Get values to insert.
    j = 0;
    for (const MDObject* obj : row) {
        addLabel(obj->label);
        mdValues[j] = row.getObject(obj->label);
        j++;
    }
    mdValues.resize(j);

    // Execute statement.
    if (!myMDSql->setObjectValues( -1, mdValues))
    {
        std::cerr << "execAddRow: error executing myMDSql->setObjectValues" << std::endl;
        success = false;
    }

    return(success);
}

void MetaDataDb::finalizeAddRow(void)
{
    myMDSql->finalizePreparedStmt();
}

size_t MetaDataDb::addRow(const MDRow &row)
{
    size_t id = addObject();
    for (auto obj : row)
        if (obj->label != MDL_FIRST_LABEL)
            setValue(*obj, id);

    return id;
}

bool MetaDataDb::getRowValues(size_t id, std::vector<MDObject> &values) const {
    for (auto &v : values) {
        if (!containsLabel(v.label))
                return false;
    }
    if (id == BAD_OBJID)
        REPORT_ERROR(ERR_MD_NOACTIVE, "getValue: please provide objId other than -1");
    return sqlUtils::select(id,
            myMDSql->db,
            myMDSql->tableName(myMDSql->tableId),
            values);
}

void MetaDataDb::addRowOpt(const MDRowSql &row)
{
    addRows({row});
}

void MetaDataDb::addMissingLabels(const MDRow &row) {
    // find missing labels
    std::vector<MDLabel> missingLabels;
    auto definedLabels = row.labels();
    for (const auto &l : definedLabels){
        if ( ! containsLabel(l)) {
            missingLabels.emplace_back(l);
        }
    }
    // add missing labels
    if ( ! missingLabels.empty()) {
        sqlUtils::addColumns(missingLabels,
                    myMDSql->db,
                    myMDSql->tableName(myMDSql->tableId));
        this->_activeLabels.insert(this->_activeLabels.end(), missingLabels.begin(), missingLabels.end());
    }
}

void MetaDataDb::addRows(const std::vector<MDRowSql> &rows)
{
    const auto noOfRows = rows.size();
    if (0 == noOfRows) {
        return;
    }
    const auto &firstRow = rows.at(0);

    // assuming all rows are using the same labels
    addMissingLabels(firstRow);

    // create mask of valid labels
    std::vector<MDLabel> labels;
    labels.reserve(firstRow.size());
    for (const MDObject* obj : firstRow)
        labels.emplace_back(obj->label);
    const auto noOfLabels = labels.size();

    // extract values to be added
    std::vector<std::vector<const MDObject*>> records;
    records.reserve(noOfRows);
    for (const auto &r : rows) {
        records.emplace_back(std::vector<const MDObject*>());
        auto &vals = records.back();
        vals.reserve(noOfLabels);
        for (const auto &l : labels) {
            vals.emplace_back(r.getObject(l));
        }
    }
    // insert values to db
    sqlUtils::insert(records, myMDSql->db,
                    myMDSql->tableName(myMDSql->tableId));
}


size_t MetaDataDb::addRow2(const MDRow &row)
{
    size_t id = BAD_OBJID;

    // Initialize INSERT.
    if (initAddRow( row))
    {
        // Execute INSERT.
        if (execAddRow( row))
        {
            // Get last inserted row id.
            id = myMDSql->getObjId();
        }

        // Finalize INSERT.
        finalizeAddRow();
    }

    return(id);
}

MetaDataDb::MetaDataDb()
{
    myMDSql = new MDSql(this);
    init({});
}//close MetaData default Constructor

MetaDataDb::MetaDataDb(const MetaData &md) {
    myMDSql = new MDSql(this);
    init({});
    MetaData::operator=(md);
}

MetaDataDb::MetaDataDb(const std::vector<MDLabel> &labelsVector)
{
    myMDSql = new MDSql(this);
    init(labelsVector);
}//close MetaData default Constructor

MetaDataDb::MetaDataDb(const FileName &fileName, const std::vector<MDLabel> &desiredLabels)
{
    myMDSql = new MDSql(this);
    init(desiredLabels);
    read(fileName, desiredLabels.empty() ? nullptr : &desiredLabels);
}//close MetaData from file Constructor

MetaDataDb::MetaDataDb(const MetaDataDb &md)
{
    myMDSql = new MDSql(this);
    copyMetadata(md);
}//close MetaData copy Constructor

MetaDataDb& MetaDataDb::operator=(const MetaDataDb &md)
{
    copyMetadata(md);
    return *this;
}

MetaDataDb::~MetaDataDb()
{
    _clear();
    delete myMDSql;
}//close MetaData Destructor

//-------- Getters and Setters ----------

int MetaDataDb::getMaxStringLength(const MDLabel thisLabel) const
{
    if (!containsLabel(thisLabel))
        return -1;

    return myMDSql->columnMaxLength(thisLabel);
}

size_t MetaDataDb::size() const
{
    return myMDSql->size();
}

bool MetaDataDb::addLabel(const MDLabel label, int pos)
{
    if (containsLabel(label))
        return false;
    if (pos < 0 || pos >= (int)this->_activeLabels.size())
        this->_activeLabels.emplace_back(label);
    else
        this->_activeLabels.insert(this->_activeLabels.begin() + pos, label);
    myMDSql->addColumn(label);
    return true;
}

bool MetaDataDb::removeLabel(const MDLabel label)
{
    std::vector<MDLabel>::iterator location;
    location = std::find(this->_activeLabels.begin(), this->_activeLabels.end(), label);

    if (location == this->_activeLabels.end())
        return false;

    this->_activeLabels.erase(location);
    return true;
}

size_t MetaDataDb::addObject()
{
    return (size_t)myMDSql->addRow();
}

void MetaDataDb::importObject(const MetaData &md, const size_t id, bool doClear)
{
    // Currently supports importing only from MetaDataDb
    assert(dynamic_cast<const MetaDataDb*>(&md) != nullptr);

    const MetaDataDb& mdd = dynamic_cast<const MetaDataDb&>(md);
    MDValueEQ query(MDL_OBJID, id);
    mdd.myMDSql->copyObjects(this, &query);
}

void MetaDataDb::importObjects(const MetaData &md, const std::vector<size_t> &objectsToAdd, bool doClear)
{
    const std::vector<MDLabel>& labels = md.getActiveLabels();
    init(labels);
    copyInfo(md);
    int size = objectsToAdd.size();
    for (int i = 0; i < size; i++)
        importObject(md, objectsToAdd[i]);
}

void MetaDataDb::importObjects(const MetaData &md, const MDQuery &query, bool doClear)
{
    // Currently supports importing only from MetaDataDb
    assert(dynamic_cast<const MetaDataDb*>(&md) != nullptr);

    const MetaDataDb& mdd = dynamic_cast<const MetaDataDb&>(md);
    this->_importObjectsDb(mdd, query, doClear);
}

void MetaDataDb::_importObjectsDb(const MetaDataDb &md, const MDQuery &query, bool doClear)
{
    if (doClear)
    {
        //Copy all structure and info from the other metadata
        init(md._activeLabels);
        copyInfo(md);
    }
    else
    {
        //If not clear, ensure that the have the same labels
        for (size_t i = 0; i < md._activeLabels.size(); i++)
            addLabel(md._activeLabels[i]);
    }
    md.myMDSql->copyObjects(this, &query);
}

bool MetaDataDb::removeObject(size_t id) {
    int removed = removeObjects(MDValueEQ(MDL_OBJID, id));
    return (removed > 0);
}

void MetaDataDb::removeObjects(const std::vector<size_t> &toRemove)
{
    int size = toRemove.size();
    for (int i = 0; i < size; i++)
        removeObject(toRemove[i]);
}

int MetaDataDb::removeObjects(const MDQuery &query)
{
    int removed = myMDSql->deleteObjects(&query);
    return removed;
}

int MetaDataDb::removeObjects()
{
    int removed = myMDSql->deleteObjects();
    return removed;
}

void MetaDataDb::addIndex(MDLabel label) const
{
    std::vector<MDLabel> labels(1);
    labels[0]=label;
    addIndex(labels);
}
void MetaDataDb::addIndex(const std::vector<MDLabel> &desiredLabels) const
{

    myMDSql->indexModify(desiredLabels, true);
}

void MetaDataDb::removeIndex(MDLabel label)
{
    std::vector<MDLabel> labels(1);
    labels[0]=label;
    removeIndex(labels);
}

void MetaDataDb::removeIndex(const std::vector<MDLabel> &desiredLabels)
{
    myMDSql->indexModify(desiredLabels, false);
}

void MetaDataDb::addItemId()
{
    addLabel(MDL_ITEM_ID);
    fillLinear(MDL_ITEM_ID,1,1);
}

void MetaDataDb::removeItemId()
{
    removeLabel(MDL_ITEM_ID);
}

//----------Iteration functions -------------------

size_t MetaDataDb::firstRowId() const
{
    return myMDSql->firstRow();
}

size_t MetaDataDb::firstObject(const MDQuery & query) const
{
    std::vector<size_t> ids;
    findObjects(ids, query);
    size_t id = ids.size() == 1 ? ids[0] : BAD_OBJID;
    return id;
}

size_t MetaDataDb::lastRowId() const
{
    return myMDSql->lastRow();
}

//-------------Search functions-------------------
void MetaDataDb::findObjects(std::vector<size_t> &objectsOut, const MDQuery &query) const
{
    objectsOut.clear();
    myMDSql->selectObjects(objectsOut, &query);
}

void MetaDataDb::findObjects(std::vector<size_t> &objectsOut, int limit) const
{
    objectsOut.clear();
    MDQuery query(limit);
    myMDSql->selectObjects(objectsOut, &query);
}

size_t MetaDataDb::countObjects(const MDQuery &query) const
{
    std::vector<size_t> objects;
    findObjects(objects, query);
    return objects.size();
}

bool MetaDataDb::containsObject(size_t objectId) const
{
    return containsObject(MDValueEQ(MDL_OBJID, objectId));
}

bool MetaDataDb::containsObject(const MDQuery &query) const
{
    std::vector<size_t> objects;
    findObjects(objects, query);
    return objects.size() > 0;
}

//--------------IO functions -----------------------
#include <sys/stat.h>
#include <fcntl.h>
#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif

void MetaDataDb::write(const FileName &_outFile, WriteModeMetaData mode) const
{
    String blockName;
    FileName outFile;
    FileName extFile;

    blockName=_outFile.getBlockName();
    if (blockName.empty())
        blockName = DEFAULT_BLOCK_NAME;
    outFile = _outFile.removeBlockName();
    extFile = _outFile.getExtension();

    if (extFile=="xml")
    {
        writeXML(outFile, blockName, mode);
    }
    else if(extFile=="sqlite")
    {
        writeDB(outFile, blockName, mode);
    }
    else
    {
        writeStar(outFile, blockName, mode);
    }
}

void MetaDataDb::_writeRows(std::ostream &os) const
{

    auto sortedLabels = this->_activeLabels;
    std::sort(sortedLabels.begin(), sortedLabels.end());    
    for (const auto& row : *this)
    {
        for (size_t i = 0; i < sortedLabels.size(); i++)
        {
            if (sortedLabels[i] != MDL_STAR_COMMENT)
            {
                os.width(1);
                row.getObject(sortedLabels[i])->toStream(os, true);
                os << " ";
            }
        }

        os << '\n';
    }
}

void MetaDataDb::write(std::ostream &os,const String &blockName, WriteModeMetaData mode ) const
{
    if(mode==MD_OVERWRITE)
        os << FileNameVersion << " * "// << (isColumnFormat ? "column" : "row")
        << '\n' //write which type of format (column or row) and the path;
        << WordWrap(this->_comment, line_max);     //write md comment in the 2nd comment line of header
    //write data block
    String _szBlockName("data_");
    _szBlockName += blockName;

    if (this->isColumnFormat())
    {
        //write md columns in 3rd comment line of the header
        os << _szBlockName << '\n';
        os << "loop_" << '\n';
        auto sortedLabels = this->_activeLabels;
        std::sort(sortedLabels.begin(), sortedLabels.end());
        for (size_t i = 0; i < sortedLabels.size(); i++)
        {
            const auto &label = sortedLabels.at(i);
            if (label != MDL_STAR_COMMENT)
            {
                os << " _" << MDL::label2Str(label) << '\n';
            }
        }
        _writeRows(os);

        //Put the activeObject to the first, if exists
    }
    else //rowFormat
    {
        os << _szBlockName << '\n';

        // Get first object. In this case (row format) there is a single object
        size_t id = firstRowId();

        if (id != BAD_OBJID)
        {
            auto sortedLabels = this->_activeLabels;
            std::sort(sortedLabels.begin(), sortedLabels.end());
            for (size_t i = 0; i < sortedLabels.size(); i++)
            {
                const auto &label = sortedLabels.at(i);
                if (label != MDL_STAR_COMMENT)
                {
                    MDObject mdValue(label);
                    os << " _" << MDL::label2Str(label) << " ";
                    myMDSql->getObjectValue(id, mdValue);
                    mdValue.toStream(os);
                    os << '\n';
                }
            }
        }

    }
}//write


void MetaDataDb::_parseObjects(std::istream &is, std::vector<MDObject*> &columnValues, const std::vector<MDLabel> *desiredLabels, bool firstTime)
{
    size_t i=0;             // Loop counter.
    size_t size=0;          // Column values vector size.

    // Columns loop.
    size = columnValues.size();
    for (i=0; i<size ;i++)
    {
        columnValues[i]->fromStream(is);
        if (is.fail())
        {
           String errorMsg = formatString("MetaData: Error parsing column '%s' value.", MDL::label2Str(columnValues[i]->label).c_str());
           columnValues[i]->failed = true;
           std::cerr << "WARNING: " << errorMsg << std::endl;
           //REPORT_ERROR(ERR_MD_BADLABEL, (String)"read: Error parsing data column, expecting " + MDL::label2Str(object.label));
        }
        else
        {
            if (firstTime)
            {
                // Check if current column label exists.
                if (columnValues[i]->label != MDL_UNDEFINED)
                {
                    // If there are no desired labels then add all.
                    bool reallyAdd=false;
                    if (desiredLabels==NULL)
                    {
                        reallyAdd=true;
                    }
                    else
                    {
                        // Check if current column belongs to desired labels.
                        for (size_t j=0; j<desiredLabels->size(); ++j)
                        {
                            if ((*desiredLabels)[j]==columnValues[i]->label)
                            {
                                reallyAdd=true;
                                break;
                            }
                        }
                    }

                    // Add label if not exists.
                    if (reallyAdd)
                    {
                        addLabel(columnValues[i]->label);
                    }
                }
            }
        }
    }

    // Insert elements in DB.
    myMDSql->setObjectValues( -1, columnValues, desiredLabels);
}



void MetaDataDb::read(const FileName &_filename,
                    const std::vector<MDLabel> *desiredLabels,
                    bool decomposeStack)
{
    String blockName;
    FileName inFile;

    blockName=_filename.getBlockName();
    //    if (blockName.empty())
    //        blockName = DEFAULT_BLOCK_NAME;
    inFile = _filename.removeBlockName();
    String extFile = _filename.getExtension();
    blockName=escapeForRegularExpressions(blockName);

    _clear();
    myMDSql->createMd();
    this->setColumnFormat(true);

    if (extFile=="xml")
        readXML(inFile, desiredLabels, blockName, decomposeStack);
    else if(extFile=="sqlite")
        readDB(inFile, desiredLabels, blockName, decomposeStack);
    else
        readStar(_filename, desiredLabels, blockName, decomposeStack);
}


void MetaDataDb::readPlain(const FileName &inFile, const String &labelsString, const String &separator)
{
    constexpr size_t LINE_LENGTH = 1024;
    clear();
    std::vector<MDLabel> labels;
    MDL::str2LabelVector(labelsString, labels);

    char lineBuffer[LINE_LENGTH];
    String line;
    std::ifstream is(inFile.c_str(), std::ios_base::in);
    size_t lineCounter = 0;
    size_t columnsNumber = labels.size();
    size_t objId;
    StringVector parts;

    while (is.getline(lineBuffer, LINE_LENGTH))
    {
        ++lineCounter;
        line.assign(lineBuffer);
        trim(line);
        if (line[0]=='#') // This is an old Xmipp comment
            continue;
        if (!line.empty())
        {
            std::stringstream ss(line);
            objId = addObject();
            for (size_t i = 0; i < columnsNumber; ++i)
            {
                MDObject obj(labels[i]);
                _parseObject(ss, obj, objId);
                setValue(obj, objId);
            }
        }
    }
}

void MetaDataDb::addPlain(const FileName &inFile, const String &labelsString, const String &separator)
{
    MetaDataDb md2;
    md2.readPlain(inFile, labelsString);
    merge(md2);
}

bool MetaDataDb::existsBlock(const FileName &_inFile)
{
#ifdef XMIPP_MMAP
    String blockName;
    FileName outFile;

    blockName=_inFile.getBlockName();
    outFile = _inFile.removeBlockName();

    struct stat file_status;
    int fd;
    char *map;

    //check if file exists or not block name has been given
    //in our format no two identical data_xxx strings may exists

    if (blockName.empty() || !outFile.exists())
        return false;
    else
    {
        //does blockname exists?
        //remove it from file in this case
        // get length of file:
        if(stat(outFile.c_str(), &file_status) != 0)
            REPORT_ERROR(ERR_IO_NOPATH,"Metadata:existsBlock can not get filesize for file "+outFile);
        size_t size = file_status.st_size;
        if(size!=0)//size=0 for /dev/stderr
        {
            fd = open(outFile.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
            if (fd == -1)
                REPORT_ERROR(ERR_IO_NOPATH,"Metadata:existsBlock can not read file named "+outFile);

            map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (map == MAP_FAILED)
                REPORT_ERROR(ERR_MEM_BADREQUEST,"Metadata:existsBlock can not map memory ");

            // Is this a START formatted FILE
            String _szBlockName = (String)("\ndata_") + blockName;
            size_t blockNameSize = _szBlockName.size();
            close(fd);
            bool found=_memmem(map, size, _szBlockName.data(), blockNameSize) != NULL;
            if (munmap(map, size) == -1)
                REPORT_ERROR(ERR_MEM_NOTDEALLOC,"metadata:write, Can not unmap memory");
            return found;
        }
        return false;
    }
#else
    REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif
}
void MetaDataDb::readXML(const FileName &filename,
                       const std::vector<MDLabel> *desiredLabels,
                       const String & blockRegExp,
                       bool decomposeStack)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"readXML not implemented yet");
}

void MetaDataDb::readDB(const FileName &filename,
                      const std::vector<MDLabel> *desiredLabels,
                      const String & blockRegExp,
                      bool decomposeStack)//what is decompose stack for?
{
    myMDSql->copyTableFromFileDB(blockRegExp, filename, desiredLabels, _maxRows);
}

/* This function parses rows data in START format
 */
void MetaDataDb::_readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues,
                               const std::vector<MDLabel> *desiredLabels) {
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

    if (myMDSql->initializeInsert( desiredLabels, columnValues))
    {
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

        myMDSql->finalizePreparedStmt();
    }

    delete[] buffer;
}

/*This function will read the md data if is in row format */
void MetaDataDb::_readRowFormat(std::istream& is) {
    String line, token;
    MDLabel label;

    size_t objectID = addObject();

    // Read data and fill structures accordingly
    while (getline(is, line, '\n'))
    {
        if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
            continue;

        // Parse labels
        std::stringstream os(line);

        os >> token;
        label = MDL::str2Label(token);
        MDObject value(label);
        os >> value;
        if (label != MDL_UNDEFINED)
            setValue(value, objectID);
    }
}

void MetaDataDb::merge(const MetaData &md2)
{
    if (size() != md2.size())
        REPORT_ERROR(ERR_MD, "Size of two metadatas should coincide for merging.");

    for (const auto& row : md2)
        this->setRow(row, row.id());
}

#define SET_AND_FILL() generator.label=label; generator.fill(*this)

void MetaDataDb::fillExpand(MDLabel label)
{
    //aggregate metadata by label (that is, avoid repetitions
    MetaDataDb mdCTFs;
    mdCTFs.distinct(*this,label);
    //read file-metadatas in new metadata
    MetaDataDb ctfModel;
    FileName fn;
    MDRowSql row;

    for (size_t id : mdCTFs.ids())
    {
        if (mdCTFs.getValue(label, fn, id))
        {
            ctfModel.read(fn);
            if (ctfModel.isEmpty())
                REPORT_ERROR(ERR_VALUE_INCORRECT, "Only can expand non empty metadatas");
            ctfModel.getRow(row, ctfModel.firstRowId());
            mdCTFs.setRow(row, id);
        }
    }
    //join
    MetaDataDb md(*this);
    join1(md, mdCTFs, label);
}

void MetaDataDb::fillConstant(MDLabel label, const String &value)
{
    MDConstGenerator generator(value);
    SET_AND_FILL();
}

void MetaDataDb::fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3)
{
    MDRandGenerator generator(op1, op2, mode, op3);
    SET_AND_FILL();
}

void MetaDataDb::fillLinear(MDLabel label, double initial, double step)
{
    MDLinealGenerator generator(initial, step);
    SET_AND_FILL();
}

void MetaDataDb::copyColumn(MDLabel labelDest, MDLabel labelSrc)
{
    String srcName = MDL::label2Str(labelSrc);
    if (!containsLabel(labelSrc))
        REPORT_ERROR(ERR_ARG_MISSING, formatString("Source label: '%s' doesn't exist on metadata", srcName.c_str()));
    addLabel(labelDest);

    String destName = MDL::label2Str(labelDest);
    String cmd = formatString("%s=%s", destName.c_str(), srcName.c_str());
    operate(cmd);
}

void MetaDataDb::copyColumnTo(MetaData &md, MDLabel labelDest, MDLabel labelSrc)
{
    if (!containsLabel(labelSrc))
        REPORT_ERROR(ERR_ARG_MISSING, formatString("Source label: '%s' doesn't exist on metadata",
                     (MDL::label2Str(labelSrc)).c_str()));
    md.addLabel(labelDest);
    std::vector<MDObject> values;
    getColumnValues(labelSrc, values);
    md.setColumnValues(values);
}

void MetaDataDb::renameColumn(MDLabel oldLabel, MDLabel newLabel)
{
    if (!containsLabel(oldLabel))
        REPORT_ERROR(ERR_ARG_MISSING, formatString("Source label: '%s' doesn't exist on metadata",
                     (MDL::label2Str(oldLabel)).c_str()));
    std::vector<MDLabel> vOldLabel(1);
    vOldLabel[0]=oldLabel;
    std::vector<MDLabel> vNewLabel(1);
    vNewLabel[0]=newLabel;
    renameColumn(vOldLabel,vNewLabel);
}

void MetaDataDb::renameColumn(const std::vector<MDLabel> &vOldLabel,
                            const std::vector<MDLabel> &vNewLabel)
{
    myMDSql->renameColumn(vOldLabel,vNewLabel);
}

void MetaDataDb::aggregateSingle(MDObject &mdValueOut, AggregateOperation op,
                               MDLabel aggregateLabel)

{
    mdValueOut.setValue(myMDSql->aggregateSingleDouble(op,aggregateLabel));
}

void MetaDataDb::aggregateSingleSizeT(MDObject &mdValueOut, AggregateOperation op,
                                    MDLabel aggregateLabel)

{
    mdValueOut.setValue(myMDSql->aggregateSingleSizeT(op,aggregateLabel));
}


double MetaDataDb::getColumnMax(MDLabel column)
{
    MDObject result(column);
    aggregateSingle(result, AGGR_MAX, column);
    return result.getValue2(double());
}

double MetaDataDb::getColumnMin(MDLabel column)
{
    MDObject result(column);
    aggregateSingle(result, AGGR_MIN, column);
    return result.getValue2(double());
}


void MetaDataDb::aggregateSingleInt(MDObject &mdValueOut, AggregateOperation op,
                                  MDLabel aggregateLabel)

{
    size_t aux = myMDSql->aggregateSingleSizeT(op,aggregateLabel);
    int aux2 = (int) aux;
    mdValueOut.setValue(aux2);
}

void MetaDataDb::aggregate(const MetaDataDb &mdIn, AggregateOperation op,
                         MDLabel aggregateLabel, MDLabel operateLabel, MDLabel resultLabel)
{
    std::vector<MDLabel> labels(2);
    std::vector<MDLabel> operateLabels(1);
    labels[0] = aggregateLabel;
    labels[1] = resultLabel;
    operateLabels[0]=operateLabel;
    init(labels);
    std::vector<AggregateOperation> ops(1);
    ops[0] = op;
    mdIn.myMDSql->aggregateMd(this, ops, operateLabels);
}

void MetaDataDb::aggregate(const MetaDataDb &mdIn, const std::vector<AggregateOperation> &ops,
                         const std::vector<MDLabel> &operateLabels,
                         const std::vector<MDLabel> &resultLabels)
{
    if (resultLabels.size() - ops.size() != 1)
        REPORT_ERROR(ERR_MD, "Labels vectors should contain one element more than operations");
    init(resultLabels);
    mdIn.myMDSql->aggregateMd(this, ops, operateLabels);
}

void MetaDataDb::aggregateGroupBy(const MetaDataDb &mdIn,
                                AggregateOperation op,
                                const std::vector<MDLabel> &groupByLabels,
                                MDLabel operateLabel,
                                MDLabel resultLabel)
{
    std::vector<MDLabel> labels;
    labels = groupByLabels;
    labels.emplace_back(resultLabel);
    init(labels);
    mdIn.myMDSql->aggregateMdGroupBy(this, op, groupByLabels, operateLabel, resultLabel);
}

//-------------Set Operations ----------------------
void MetaDataDb::_setOperates(const MetaDataDb &mdIn,
                            const MDLabel label,
                            SetOperation operation)
{
    std::vector<MDLabel> labels;
    labels.emplace_back(label);
    _setOperates(mdIn,labels,operation);
}

void MetaDataDb::_setOperates(const MetaDataDb &mdIn,
                            const std::vector<MDLabel> &labels,
                            SetOperation operation)
{
    if (this == &mdIn) //not sense to operate on same metadata
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation on input metadata");
    if (size() == 0 && mdIn.size() == 0)
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation if both metadata are empty");
    //Add labels to be sure are present
    for (size_t i = 0; i < mdIn._activeLabels.size(); i++)
        addLabel(mdIn._activeLabels[i]);

    mdIn.myMDSql->setOperate(this, labels, operation);
}

void MetaDataDb::_setOperatesLabel(const MetaDataDb &mdIn,
                            const MDLabel label,
                            SetOperation operation)
{
    if (this == &mdIn) //not sense to operate on same metadata
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation on input metadata");
    if (mdIn.size() == 0)
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation if both metadata are empty");
    //Add label to be sure is present in output
    addLabel(label);
    std::vector<MDLabel> labels;
    labels.emplace_back(label);
    mdIn.myMDSql->setOperate(this, labels, operation);
}

void MetaDataDb::_setOperates(const MetaDataDb &mdInLeft,
                            const MetaDataDb &mdInRight,
                            const std::vector<MDLabel> &labelsLeft,
                            const std::vector<MDLabel> &labelsRight,
                            SetOperation operation)
{
    if (this == &mdInLeft || this == &mdInRight) //not sense to operate on same metadata
        REPORT_ERROR(ERR_MD, "Couldn't perform this operation on input metadata");
    //Add labels to be sure are present
    for (size_t i = 0; i < mdInLeft._activeLabels.size(); i++)
        addLabel(mdInLeft._activeLabels[i]);
    for (size_t i = 0; i < mdInRight._activeLabels.size(); i++)
    {
        bool found=false;
        for (size_t j=0; j<labelsRight.size(); ++j)
            if (mdInRight._activeLabels[i]==labelsRight[j])
            {
                found=true;
                break;
            }
        if (!found)
            addLabel(mdInRight._activeLabels[i]);
    }

    myMDSql->setOperate(&mdInLeft, &mdInRight, labelsLeft,labelsRight, operation);
}

void MetaDataDb::unionDistinct(const MetaDataDb &mdIn, const MDLabel label)
{
    if(mdIn.isEmpty())
        return;
    _setOperates(mdIn, label, UNION_DISTINCT);
}

void MetaDataDb::unionAll(const MetaDataDb &mdIn)
{
    if(mdIn.isEmpty())
        return;
    _setOperates(mdIn, MDL_UNDEFINED, UNION);//label not needed for unionAll operation
}


void MetaDataDb::intersection(const MetaDataDb &mdIn, const MDLabel label)
{
    if(mdIn.isEmpty())
        clear();
    else
        _setOperates(mdIn, label, INTERSECTION);
}

void MetaDataDb::removeDuplicates(MetaDataDb &MDin, MDLabel label)
{
    if(MDin.isEmpty())
        return;
    _setOperates(MDin, label, REMOVE_DUPLICATE);
}

void MetaDataDb::distinct(MetaDataDb &MDin, MDLabel label)
{
    if(MDin.isEmpty())
        return;
    _setOperatesLabel(MDin, label, DISTINCT);
}

void MetaDataDb::subtraction(const MetaDataDb &mdIn, const MDLabel label)
{
    if(mdIn.isEmpty())
        return;
    _setOperates(mdIn, label, SUBSTRACTION);
}

void MetaDataDb::join1(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const MDLabel label, JoinType type)
{
    join2(mdInLeft, mdInRight, label, label, type);
}

void MetaDataDb::join2(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const MDLabel labelLeft,
                    const MDLabel labelRight, JoinType type)
{
    clear();
    std::vector<MDLabel> labelsLeft, labelsRight;
    labelsLeft.emplace_back(labelLeft);
    labelsRight.emplace_back(labelRight);
    _setOperates(mdInLeft, mdInRight, labelsLeft,labelsRight, (SetOperation)type);
}

void MetaDataDb::join1(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const std::vector<MDLabel> &labels, JoinType type)
{
    join2(mdInLeft, mdInRight, labels, labels, type);
}

void MetaDataDb::join2(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const std::vector<MDLabel> &labelsLeft,
                    const std::vector<MDLabel> &labelsRight, JoinType type)
{
    clear();
    _setOperates(mdInLeft, mdInRight, labelsLeft,labelsRight, (SetOperation)type);
}

void MetaDataDb::joinNatural(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight)
{
    join2(mdInLeft, mdInRight, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL);
}

void MetaDataDb::operate(const String &expression)
{
    if (!myMDSql->operate(expression))
        REPORT_ERROR(ERR_MD, "MetaDataDb::operate: error doing operation");
}

void MetaDataDb::replace(const MDLabel label, const String &oldStr, const String &newStr)
{
    String labelStr = MDL::label2Str(label);
    String expression = formatString("%s=replace(%s,'%s', '%s')",
                                     labelStr.c_str(), labelStr.c_str(), oldStr.c_str(), newStr.c_str());
    if (!myMDSql->operate(expression))
        REPORT_ERROR(ERR_MD, "MetaDataDb::replace: error doing operation");
}

void MetaDataDb::randomize(const MetaDataDb &MDin)
{
    std::random_device rd;
    auto g = std::mt19937(rd());
    std::vector<size_t> objects;
    MDin.myMDSql->selectObjects(objects);
    std::shuffle(objects.begin(), objects.end(), g);
    importObjects(MDin, objects);
}

void MetaDataDb::sort(MetaDataDb &MDin, const MDLabel sortLabel,bool asc, int limit, int offset)
{
    if (MDin.containsLabel(sortLabel))
    {
        init(MDin._activeLabels);
        copyInfo(MDin);
        //if you sort just once the index will not help much
        addIndex(sortLabel);
        MDQuery query(limit, offset, sortLabel,asc);
        MDin.myMDSql->copyObjects(this, &query);
    }
    else
        *this=MDin;
}

void MetaDataDb::sort(MetaDataDb &MDin, const String &sortLabel,bool asc, int limit, int offset)
{
    // Check if the label has semicolon
    size_t ipos=sortLabel.find(':');
    MDLabelType type = MDL::labelType(sortLabel);
    if (ipos!=String::npos || type == LABEL_VECTOR_DOUBLE || type == LABEL_VECTOR_SIZET)
    {
        if(limit != -1 || offset != 0)
            REPORT_ERROR(ERR_ARG_INCORRECT,"Limit and Offset are not implemented for vector sorting.");

        MDLabel label;
        size_t column;
        if (ipos!=String::npos)
        {
            // Check that the label is a vector field
            std::vector< String > results;
            splitString(sortLabel,":",results);
            column=textToInteger(results[1]);
            MDLabelType type = MDL::labelType(results[0]);
            if (type != LABEL_VECTOR_DOUBLE || type != LABEL_VECTOR_SIZET)
                REPORT_ERROR(ERR_ARG_INCORRECT,"Column specifications cannot be used with non-vector labels");
            label = MDL::str2Label(results[0]);
        }
        else
        {
            label = MDL::str2Label(sortLabel);
            column = 0;
        }

        // Get the column values
        MultidimArray<double> v;
        v.resizeNoCopy(MDin.size());
        std::vector<double> vectorValues;
        int i = 0;
        for (size_t id : MDin.ids())
        {
            MDin.getValue(label, vectorValues, id);
            if (column >= vectorValues.size())
                REPORT_ERROR(ERR_MULTIDIM_SIZE,"Trying to access to inexistent column in vector");
            DIRECT_A1D_ELEM(v, i) = vectorValues[column];
            i++;
        }

        // Sort
        MultidimArray<int> idx;
        v.indexSort(idx);

        // Construct output Metadata
        init(MDin._activeLabels);
        copyInfo(MDin);
        size_t id;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(idx)
        {
            MDRowSql row;
            MDin.getRow(row,DIRECT_A1D_ELEM(idx,i));
            id = addObject();
            setRow(row, id);
        }
    }
    else
    {
        sort(MDin, MDL::str2Label(sortLabel),asc, limit, offset);
    }
}

void MetaDataDb::split(size_t n, std::vector<MetaDataDb> &results, const MDLabel sortLabel)
{
    size_t mdSize = size();
    if (n > mdSize)
        REPORT_ERROR(ERR_MD, "MetaDataDb::split: Couldn't split a metadata in more parts than its size");

    results.clear();
    results.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        MetaDataDb &md = results.at(i);
        md._selectSplitPart(*this, n, i, mdSize, sortLabel);
    }
}

void MetaDataDb::_selectSplitPart(const MetaDataDb &mdIn,
                                int n, int part, size_t mdSize,
                                const MDLabel sortLabel)
{
    size_t first, last, n_images;
    n_images = divide_equally(mdSize, n, part, first, last);
    init(mdIn._activeLabels);
    copyInfo(mdIn);
    mdIn.myMDSql->copyObjects(this, new MDQuery(n_images, first, sortLabel));
}

void MetaDataDb::selectSplitPart(const MetaData &mdIn, size_t n, size_t part, const MDLabel sortLabel)
{
    if (dynamic_cast<const MetaDataDb*>(&mdIn) != nullptr)
        return _selectSplitPart(dynamic_cast<const MetaDataDb&>(mdIn), n, part, sortLabel);
    throw std::logic_error("Not yet implemented"); // TODO: use universal functions just on MetaData
}

void MetaDataDb::_selectSplitPart(const MetaDataDb &mdIn, size_t n, size_t part, const MDLabel sortLabel)
{
    size_t mdSize = mdIn.size();
    if (n > mdSize)
        REPORT_ERROR(ERR_MD, "selectSplitPart: Couldn't split a metadata in more parts than its size");
    if (part < 0 || part >= n)
        REPORT_ERROR(ERR_MD, "selectSplitPart: 'part' should be between 0 and n-1");
    _selectSplitPart(mdIn, n, part, mdSize, sortLabel);

}

void MetaDataDb::selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel)
{
    if (dynamic_cast<const MetaDataDb*>(&mdIn) != nullptr)
        return _selectRandomSubset(dynamic_cast<const MetaDataDb&>(mdIn), numberOfObjects, sortLabel);
    throw std::logic_error("Not yet implemented"); // TODO: use universal functions just on MetaData
}

void MetaDataDb::_selectRandomSubset(const MetaDataDb &mdIn, size_t numberOfObjects, const MDLabel sortLabel)
{
    clear();

    MetaDataDb mdAux, mdAux2;
    mdAux.randomize(mdIn);
    mdAux2.selectPart(mdAux, 0, numberOfObjects);
    sort(mdAux2,sortLabel);
}

void MetaDataDb::selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                            const MDLabel sortLabel)
{
    if (dynamic_cast<const MetaDataDb*>(&mdIn) != nullptr)
        return _selectPart(dynamic_cast<const MetaDataDb&>(mdIn), startPosition, numberOfObjects, sortLabel);
    throw std::logic_error("Not yet implemented"); // TODO: use universal functions just on MetaData
}

void MetaDataDb::_selectPart(const MetaDataDb &mdIn, size_t startPosition, size_t numberOfObjects,
                            const MDLabel sortLabel)
{
    size_t mdSize = mdIn.size();
    if (startPosition < 0 || startPosition >= mdSize)
        REPORT_ERROR(ERR_MD, "selectPart: 'startPosition' should be between 0 and size()-1");
    init(mdIn._activeLabels);
    copyInfo(mdIn);
    mdIn.myMDSql->copyObjects(this, new MDQuery(numberOfObjects, startPosition, sortLabel));
}

void MetaDataDb::makeAbsPath(const MDLabel label)
{

    String aux_string;
    String aux_string_path;
    char buffer[1024];

    if (!getcwd(buffer, 1023))
        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot get the current directory");
    String path_str(buffer);
    path_str += "/";
    getValue(label, aux_string, firstRowId());

    if (aux_string[0] == '/')
        return;

    FileName auxFile;
    for (size_t id : this->ids())
    {
        aux_string_path = path_str;
        getValue(label, auxFile, id);

        if (auxFile.isInStack())
        {
            size_t id = auxFile.find('@',0);
            auxFile.insert(id+1,aux_string_path);
            setValue(label, auxFile, id);
        }
        else
        {
            auxFile.addPrefix(aux_string_path);
            setValue(label, auxFile, id);
        }
    }
}

void MetaDataDb::writeDB(const FileName fn, const FileName blockname, WriteModeMetaData mode) const
{
    if(mode==MD_OVERWRITE)
        unlink(fn.c_str());
    myMDSql->copyTableToFileDB(blockname,fn);
}

void MetaDataDb::writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const
{
    //fixme
    ////THIS SHOULD BE IMPLEMENTED USING AN XML LIBRARY THAT HANDLES THE FILE PROPERLY
    if(mode!=MD_OVERWRITE)
        REPORT_ERROR(ERR_NOT_IMPLEMENTED,"XML is only implemented for overwrite mode");
    std::ofstream ofs(fn.c_str(), std::ios_base::out|std::ios_base::trunc);
    size_t size = this->_activeLabels.size();
    ofs <<  "<" << blockname << ">"<< '\n';
    for (size_t id : this->ids())
    {
        ofs <<  "<ROW ";
        for (size_t i = 0; i < size; i++)
        {
            if (this->_activeLabels[i] != MDL_STAR_COMMENT)
            {
                ofs << MDL::label2Str(this->_activeLabels[i]) << "=\"";
                MDObject mdValue(this->_activeLabels[i]);
                //ofs.width(1);
                myMDSql->getObjectValue(id, mdValue);
                mdValue.toStream(ofs, true);
                ofs << "\" ";
            }
        }
        ofs <<  " />" << '\n';
    }
    ofs <<  "</" << blockname << ">"<< '\n';
}

void MetaDataDb::writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const
{
    std::ofstream ofs(fn.c_str(), std::ios_base::trunc|std::ios_base::out);

    if (desiredLabels != NULL)
    {
        MetaDataDb mdAux(*this);
        mdAux._activeLabels = *desiredLabels;
        mdAux._writeRows(ofs);
    }
    else
        _writeRows(ofs);
    ofs.close();
}

void MetaDataDb::metadataToVec(std::vector<MDRowSql> &vd)
{
    for (const auto& row : *this)
        vd.emplace_back(dynamic_cast<const MDRowSql&>(row));
}

void MetaDataDb::vecToMetadata(const std::vector<MDRow> &rowMetadata)
{
    const MDRowSql row;

    for (size_t i=0;i<rowMetadata.size();i++)
        this->addRow(rowMetadata[i]);
}

bool MetaDataDb::operator==(const MetaDataDb& op) const
{
    return myMDSql->equals(*(op.myMDSql));
}

std::ostream& operator<<(std::ostream& o, const MetaData & mD)
{
    mD.write(o);
    return o;
}
