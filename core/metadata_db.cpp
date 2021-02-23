/**************************************************************************
 *
 * Authors:      J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

// Get the blocks available
void getBlocksInMetaDataFile(const FileName &inFile, StringVector& blockList)
{
    if (!inFile.isMetaData())
        return;
    if (inFile.getBlockName()!="")
        return;
    blockList.clear();
    String extFile = inFile.getExtension();

    if (extFile=="xml")
        REPORT_ERROR(ERR_NOT_IMPLEMENTED,"getBlocksInMetaDataFile");
    else if(extFile=="sqlite")
    {
        getBlocksInMetaDataFileDB(inFile,blockList);
    }
    else
    {    //map file
        int fd;
        MetaDataDb mdAux;
        mdAux.setMaxRows(1);
        mdAux.read(inFile);
        BUFFER_CREATE(bufferMap);
        mapFile(inFile, bufferMap.begin, bufferMap.size, fd);
        BUFFER_COPY(bufferMap, buffer);
        BLOCK_CREATE(block);
        String blockName;
        while (mdAux.nextBlock(buffer, block))
        {
            BLOCK_NAME(block, blockName);
            blockList.push_back(blockName);
        }

        unmapFile(bufferMap.begin, bufferMap.size, fd);
    }
}

// Does the blocks exist
bool existsBlockInMetaDataFile(const FileName &inFileWithBlock)
{
    return existsBlockInMetaDataFile(inFileWithBlock.removeBlockName(),
                                     inFileWithBlock.getBlockName());
}
bool existsBlockInMetaDataFile(const FileName &inFile, const String& inBlock)
{
    if (!inFile.isMetaData())
        return false;
    if (!inFile.getBlockName().empty())
        return inBlock == inFile.getBlockName();

    MetaDataDb MDaux(inFile);
    //map file
    int fd;
    BUFFER_CREATE(bufferMap);
    mapFile(inFile, bufferMap.begin, bufferMap.size, fd);
    BUFFER_COPY(bufferMap, buffer);
    BLOCK_CREATE(block);
    String blockName;
    bool result = false;
    while (MDaux.nextBlock(buffer, block))
    {
        BLOCK_NAME(block, blockName);
        if (inBlock == blockName)
        {
            result = true;
            break;
        }
    }

    unmapFile(bufferMap.begin, bufferMap.size, fd);
    return result;

}

//-----Constructors and related functions ------------
void MetaDataDb::_clear(bool onlyData)
{
    if (onlyData)
    {
        myMDSql->deleteObjects();
    }
    else
    {
        _path.clear();
        _comment.clear();
        _fastStringSearch.clear();
        _fastStringSearchLabel = MDL_UNDEFINED;

        _activeLabels.clear();
        _ignoreLabels.clear();
        _isColumnFormat = true;
        _inFile = FileName();
        myMDSql->clearMd();
    }
    eFilename="";
}//close clear

void MetaDataDb::clear()
{
    //_clear(true);
    init();
}

void MetaDataDb::init(const std::vector<MDLabel> *labelsVector)
{
    _clear();
    _maxRows = 0; //by default read all rows
    _parsedLines = 0; //no parsed line;
    if (labelsVector != NULL)
        _activeLabels = *labelsVector;
    //Create table in database
    myMDSql->createMd();
    _precision = 100;
    isMetadataFile = false;
}//close init

void MetaDataDb::copyMetadata(const MetaDataDb &md, bool copyObjects)
{
    if (this == &md) //not sense to copy same metadata
        return;
    init(&(md._activeLabels));
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

void MetaDataDb::asVMetaData(VMetaData &vmdOut)
{
    vmdOut.clear();
    for (const MDRow& row : *this)
        vmdOut.push_back(dynamic_cast<const MDRowSql&>(row));
}

void MetaDataDb::fromVMetaData(VMetaData &vmdIn)
{
    clear();
    for (auto& row: vmdIn)
        addRow2(row);
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

template<typename T>
bool MetaDataDb::getColumnValuesOpt(const MDLabel label, std::vector<T> &values) const {
    if (!containsLabel(label))
            return false;
    return sqlUtils::select(label,
            myMDSql->db,
            myMDSql->tableName(myMDSql->tableId),
            values);
}

/**
 *  XXX HACK Because of the cyclic dependency between MetaData/MetaData label and MetaData SQL,
 *  this cannot be in header. So we need to explicitly instantiate it
 */
template bool MetaDataDb::getColumnValuesOpt<float>(MDLabel, std::vector<float, std::allocator<float> >&) const;
template bool MetaDataDb::getColumnValuesOpt<FileName>(MDLabel, std::vector<FileName, std::allocator<FileName> >&) const;
template bool MetaDataDb::getColumnValuesOpt<int>(MDLabel, std::vector<int, std::allocator<int> >&) const;

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

bool MetaDataDb::bindValue( size_t id) const
{
    bool success=true;

    // Prepare statement.
    if (!myMDSql->bindStatement( id))
    {
        success = false;
    }

    return(success);
}

bool MetaDataDb::initGetRow( bool addWhereClause) const
{
    bool success=true;

    // Prepare statement.
    if (!myMDSql->initializeSelect( addWhereClause, this->_activeLabels))
    {
        success = false;
    }

    return(success);
}

bool MetaDataDb::execGetRow(MDRow &row) const
{
    std::vector<MDObject> mdValues;     // Vector to store values.
    mdValues.reserve(this->_activeLabels.size());

    // Clear row.
    row.clear();

    // Execute statement.
    bool success = myMDSql->getObjectsValues(this->_activeLabels, mdValues);
    if (success) {
        // Set values in row.
        for (const auto &obj : mdValues) {
            row.setValue(obj);
        }
    }

    return(success);
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

bool MetaDataDb::getRow(MDRow &row, size_t id) const
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
    for (auto &v : values) {
        row.setValue(v);
    }
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

//TODO: could be improve in a query for update the entire row
#define SET_ROW_VALUES(row) \
    for (MDObject* obj : row) {\
        setValue(obj->label, id);\
    }\

bool MetaDataDb::initSetRow(const MDRow &row)
{
    int     j=0;                    // Loop counter.
    bool    success=true;               // Return value.
    std::vector<MDLabel>   labels;      // Columns labels.

    // Set label vector size.
    labels.resize(row.size());

    // Build labels vector:
    j = 0;
    for (const MDObject* obj : row) {
        addLabel(obj->label);
        labels[j] = obj->label;
        j++;
    }
    labels.resize(j);

    // Prepare statement.
    if (!myMDSql->initializeUpdate(labels))
    {
        success = false;
    }

    return(success);
}


bool MetaDataDb::execSetRow(const MDRow &row, size_t id)
{
    int i = 0, j = 0;
    bool success = true;
    std::vector<MDObject*> mdValues;

    // Set values vector size.
    mdValues.resize(row.size());

    // Build values vector.
    j = 0;
    for (MDObject* obj : row) {
        addLabel(obj->label);
        mdValues[i] = obj;
        j++;
    }
    mdValues.resize(j);

    // Execute statement.
    if (!myMDSql->setObjectValues(id, mdValues))
        success = false;

    return success;
}

void MetaDataDb::finalizeSetRow(void)
{
    myMDSql->finalizePreparedStmt();
}


bool MetaDataDb::setRow(const MDRow &row, size_t id)
{
    SET_ROW_VALUES(row);

    return(true);
}

bool MetaDataDb::setRow2(const MDRow &row, size_t id)
{
    bool success = true;

    // Initialize UPDATE.
    success = initSetRow( row);
    if (success)
    {
        // Execute UPDATE.
        success = execSetRow( row, id);

        // Finalize UPDATE.
        finalizeSetRow();
    }

    return(success);
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
        std::cout << "initAddRow: error executing myMDSql->initializeInsert" << std::endl;
        success = false;
    }

    return success;
}


bool MetaDataDb::execAddRow(const MDRow &row)
{
    int j = 0;
    bool success = true;
    std::vector<MDObject*> mdValues;

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
        std::cout << "execAddRow: error executing myMDSql->setObjectValues" << std::endl;
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
    SET_ROW_VALUES(row);

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
            missingLabels.push_back(l);
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
        records.push_back(std::vector<const MDObject*>());
        auto &vals = records.back();
        vals.reserve(noOfLabels);
        for (const auto &l : labels) {
            vals.push_back(r.getObject(l));
        }
    }
    // insert values to db
    sqlUtils::insert(records, myMDSql->db,
                    myMDSql->tableName(myMDSql->tableId));
}


size_t MetaDataDb::addRow2(const MDRow &row)
{
    size_t id;      // Inserted row id.

    // Initialize INSERT.
    if (initAddRow( row))
    {
        // Execute INSERT.
        if (execAddRow( row))
        {
            // Get last inserted row id.
            id = myMDSql->getObjId();
        }
        else
        {
            id = BAD_OBJID;
        }

        // Finalize INSERT.
        finalizeAddRow();
    }

    return(id);
}

MetaDataDb::MetaDataDb()
{
    myMDSql = new MDSql(this);
    init(NULL);
}//close MetaData default Constructor

MetaDataDb::MetaDataDb(const std::vector<MDLabel> *labelsVector)
{
    myMDSql = new MDSql(this);
    init(labelsVector);
}//close MetaData default Constructor

MetaDataDb::MetaDataDb(const FileName &fileName, const std::vector<MDLabel> *desiredLabels)
{
    myMDSql = new MDSql(this);
    init(desiredLabels);
    read(fileName, desiredLabels);
}//close MetaData from file Constructor

MetaDataDb::MetaDataDb(const MetaDataDb &md)
{
    myMDSql = new MDSql(this);
    copyMetadata(md);
}//close MetaData copy Constructor

MetaDataDb& MetaDataDb::operator =(const MetaDataDb &md)
{
    copyMetadata(md);
    return *this;
}//close metadata operator =

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
        this->_activeLabels.push_back(label);
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

bool MetaDataDb::keepLabels(const std::vector<MDLabel> &labels)
{
    for (size_t i = 0; i < this->_activeLabels.size();)
    {
        if (!vectorContainsLabel(labels, this->_activeLabels[i]))
            removeLabel(this->_activeLabels[i]);
        else
            ++i;
    }
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
    init(&(md._activeLabels));
    copyInfo(md);
    int size = objectsToAdd.size();
    for (int i = 0; i < size; i++)
        importObject(md, objectsToAdd[i]);
}

void MetaDataDb::importObjects(const MetaDataDb &md, const MDQuery &query, bool doClear)
{
    if (doClear)
    {
        //Copy all structure and info from the other metadata
        init(&(md._activeLabels));
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

void MetaDataDb::writeStar(const FileName &outFile,const String &blockName, WriteModeMetaData mode) const
{
#ifdef XMIPP_MMAP
    if (outFile.hasImageExtension())
        REPORT_ERROR(ERR_IO,"MetaData:writeStar Trying to write metadata with image extension");

    struct stat file_status;
    int fd;
    char *map=NULL;
    char * tailMetadataFile = NULL;//auxiliary variable to keep metadata file tail in memory
    size_t size=-1;
    char * target, * target2=NULL;

    //check if file exists or not block name has been given
    //in our format no two identical data_xxx strings may exists
    if(mode == MD_APPEND)
    {
        if (blockName.empty() || !outFile.exists())
            mode = MD_OVERWRITE;
        else
        {
            //does blockname exists?
            //remove it from file in this case
            // get length of file:
            if(stat(outFile.c_str(), &file_status) != 0)
                REPORT_ERROR(ERR_IO_NOPATH,"MetaData:writeStar can not get filesize for file "+outFile);
            size = file_status.st_size;
            if (size == 0)
                mode = MD_OVERWRITE;
        }

        if (mode == MD_APPEND)//size=0 for /dev/stderr
        {
            fd = open(outFile.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
            if (fd == -1)
                REPORT_ERROR(ERR_IO_NOPATH,"MetaData:writeStar can not read file named "+outFile);

            map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (map == MAP_FAILED)
                REPORT_ERROR(ERR_MEM_BADREQUEST,"MetaData:writeStar can not map memory ");

            // Is this a metadata formatted FILE
            if(strncmp(map,FileNameVersion.c_str(),FileNameVersion.length()) !=0 )
            {
                mode = MD_OVERWRITE;
            }
            else
            {
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

    if (tailMetadataFile != NULL)
    {
        //append memory buffer to file
        //may a cat a buffer to a ofstream
        ofs.write(tailMetadataFile,(map + size) - target2);
        free(tailMetadataFile);
    }
    ofs.close();

#else

    REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

}//function writeStar

void MetaDataDb::append(const FileName &outFile) const
{
    if (outFile.exists())
    {
        std::ofstream ofs(outFile.c_str(), std::ios_base::app);
        _writeRows(ofs);
        ofs.close();
    }
    else
        write(outFile);
}

void MetaDataDb::_writeRows(std::ostream &os) const
{
    // Prepare statement.
    this->initGetRow( true);

    for (const auto& row : *this)
    {
        for (size_t i = 0; i < this->_activeLabels.size(); i++)
        {
            if (this->_activeLabels[i] != MDL_STAR_COMMENT)
            {
                os.width(1);
                row.getObject(this->_activeLabels[i])->toStream(os, true);
                os << " ";
            }
        }

        os << '\n';
    }
    // Finalize statement.
    myMDSql->finalizePreparedStmt();
}

void MetaDataDb::print() const
{
    write(std::cout);
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

    if (_isColumnFormat)
    {
        //write md columns in 3rd comment line of the header
        os << _szBlockName << '\n';
        os << "loop_" << '\n';
        const auto noOfLabels = this->_activeLabels.size();
        for (size_t i = 0; i < noOfLabels; i++)
        {
            const auto &label = this->_activeLabels.at(i);
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
            const auto noOfLabels = this->_activeLabels.size();
            for (size_t i = 0; i < noOfLabels; i++)
            {
                const auto &label = this->_activeLabels.at(i);
                if (label != MDL_STAR_COMMENT)
                {
                    MDObject mdValue(this->_activeLabels[i]);
                    os << " _" << MDL::label2Str(label) << " ";
                    myMDSql->getObjectValue(id, mdValue);
                    mdValue.toStream(os);
                    os << '\n';
                }
            }
        }

    }
}//write

/* This function will read the possible columns from the file
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels (for OLD doc files)
 */
void MetaDataDb::_readColumns(std::istream& is, std::vector<MDObject*> & columnValues,
                            const std::vector<MDLabel>* desiredLabels)
{
    String token;
    MDLabel label;

    while (is >> token)
        if (token.find('(') == String::npos)
        {
            //label is not recognized, the MDValue will be created
            //with MDL_UNDEFINED, which will be ignored while reading data
            label = MDL::str2Label(token);

            // Try to read undefined labels as String using the buffer approach
            if (label == MDL_UNDEFINED)
                label = MDL::getNewAlias(token);

            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels
            columnValues.push_back(new MDObject(label));
            if (label != MDL_UNDEFINED)
                addLabel(label);

        }
}


void MetaDataDb::_parseObjects(std::istream &is, std::vector<MDObject*> & columnValues, const std::vector<MDLabel> *desiredLabels, bool firstTime)
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


/* Helper function to parse an MDObject and set its value.
 * The parsing will be from an input stream(istream)
 * and if parsing fails, an error will be raised
 */
void MetaDataDb::_parseObject(std::istream &is, MDObject &object, size_t id)
{
    object.fromStream(is);
    if (is.fail())
    {
       String errorMsg = formatString("MetaData: Error parsing column '%s' value.", MDL::label2Str(object.label).c_str());
       object.failed = true;
       std::cerr << "WARNING: " << errorMsg << std::endl;
       //REPORT_ERROR(ERR_MD_BADLABEL, (String)"read: Error parsing data column, expecting " + MDL::label2Str(object.label));
    }
    else
        if (object.label != MDL_UNDEFINED)
            setValue(object, id);
}//end of function parseObject

//#define END_OF_LINE() ((char*) memchr (iter, '\n', end-iter+1))
#define END_OF_LINE() ((char*) memchr (iter, '\n', end-iter))

/* This function will read the possible columns from the file
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels (for new STAR files)
 */
void MetaDataDb::_readColumnsStar(mdBlock &block,
                                std::vector<MDObject*> & columnValues,
                                const std::vector<MDLabel>* desiredLabels,
                                bool addColumns,
                                size_t id)
{
    char * end = block.end;
    char * newline = NULL;
    bool found_column;
    MDLabel label;
    char * iter = block.loop;
    if (!_isColumnFormat)
    {
        iter = block.begin;
        iter = END_OF_LINE() + 1; //this should point at first label, after data_XXX
    }

    do
    {
        found_column = false;
        while (iter[0] == '#') //Skip comment
            iter = END_OF_LINE() + 1;

        //trim spaces and newlines at the beginning
        while ( isspace(iter[0]))
            ++iter;

        if (iter < end && iter[0] == '_')
        {
            found_column = true;
            ++iter; //shift _
            std::stringstream ss;
            newline = END_OF_LINE();
            //Last label and no data needs this check
            if (newline == NULL)
                newline = end;
            String s(iter, newline - iter);//get current line
            ss.str(s);//set the string of the stream
            //Take the first token which is the label
            //if the label contain spaces will fail
            ss >> s; //get the first token, the label
            label = MDL::str2Label(s);
            if (label == MDL_UNDEFINED)
                label = MDL::getNewAlias(s);

            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels

            if (label != MDL_UNDEFINED)
                addLabel(label);

            if (addColumns)
            {
                MDObject * _mdObject = new MDObject(label);
                columnValues.push_back(_mdObject);//add the value here with a char
                if(!_isColumnFormat)
                    _parseObject(ss, *_mdObject, id);
            }
            iter = newline + 1;//go to next line character
        }
    }
    while (found_column)
        ;

    // This condition fails for empty blocks
    // if (iter < block.end)
    if (iter <= block.end +1)
        block.loop = iter; //Move loop pointer to position of last found column
}

/* This function will be used to parse the rows data
 * having read the columns labels before and setting which are desired
 * the useCommentAsImage is for compatibility with old DocFile format
 * where the image were in comments
 */
void MetaDataDb::_readRows(std::istream& is, std::vector<MDObject*> & columnValues, bool useCommentAsImage)
{
    String line = "";
    while (!is.eof() && !is.fail())
    {
        //Move until the ';' or the first alphanumeric character
        while (is.peek() != ';' && isspace(is.peek()) && !is.eof())
            is.ignore(1);
        if (!is.eof())
        {
            if (is.peek() == ';')//is a comment
            {
                is.ignore(1); //ignore the ';'
                getline(is, line);
                trim(line);
            }
            else if (!isspace(is.peek()))
            {
                size_t id = addObject();
                if (line != "")//this is for old format files
                {
                    if (!useCommentAsImage)
                        setValue(MDObject(MDL_STAR_COMMENT, line), id);
                    else
                        setValue(MDObject(MDL_IMAGE, line), id);
                }
                int nCol = columnValues.size();
                for (int i = 0; i < nCol; ++i)
                    _parseObject(is, *(columnValues[i]), id);
            }
        }
    }
}

/* This function will be used to parse the rows data in START format
 */
void MetaDataDb::_readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues, const std::vector<MDLabel> *desiredLabels)
{
    String line;
    std::stringstream ss;
    size_t n = block.end - block.loop;
    bool firstTime=true;

    if (n==0)
        return;

    char * buffer = new char[n];
    memcpy(buffer, block.loop, n);
    char *iter = buffer, *end = iter + n, * newline = NULL;
    _parsedLines = 0; //Check how many lines the md have

    if (myMDSql->initializeInsert( desiredLabels, columnValues))
    {
        while (iter < end) //while there are data lines
        {
            //Assing \n position and check if NULL at the same time
            if (!(newline = END_OF_LINE()))
                newline = end;
            line.assign(iter, newline - iter);
            trim(line);

            if (!line.empty() && line[0] != '#')
            {
                //_maxRows would be > 0 if we only want to read some
                // rows from the md for performance reasons...
                // anyway the number of lines will be counted in _parsedLines
                if (_maxRows == 0 || _parsedLines < _maxRows)
                {
                    std::stringstream ss(line);
                    _parseObjects( ss, columnValues, desiredLabels, firstTime);
                    firstTime=false;
                }
                _parsedLines++;
            }
            iter = newline + 1; //go to next line
        }

        // Finalize statement.
        myMDSql->finalizePreparedStmt();
    }

    delete[] buffer;
}

/*This function will read the md data if is in row format */
void MetaDataDb::_readRowFormat(std::istream& is)
{
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
    extFile = _filename.getExtension();
    blockName=escapeForRegularExpressions(blockName);

    _clear();
    myMDSql->createMd();
    _isColumnFormat = true;

    if (extFile=="xml")
        readXML(inFile, desiredLabels, blockName, decomposeStack);
    else if(extFile=="sqlite")
        readDB(inFile, desiredLabels, blockName, decomposeStack);
    else
        readStar(_filename, desiredLabels, blockName, decomposeStack);

    //_read(filename,desiredLabels,BlockName,decomposeStack);
    //_read calls clean so I cannot use eFilename as filename ROB
    // since eFilename is reset in clean
    eFilename = _filename;
}

#define LINE_LENGTH 1024
void MetaDataDb::readPlain(const FileName &inFile, const String &labelsString, const String &separator)
{

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
void MetaDataDb::readStar(const FileName &filename,
                        const std::vector<MDLabel> *desiredLabels,
                        const String & blockRegExp,
                        bool decomposeStack)
{
    //First try to open the file as a metadata
    size_t id;
    FileName inFile = filename.removeBlockName();

    if (!(isMetadataFile = inFile.isMetaData()))//if not a metadata, try to read as image or stack
    {
        Image<char> image;
        if (decomposeStack) // If not decomposeStack it is no necessary to read the image header
            image.read(filename, HEADER);
        if ( !decomposeStack || image().ndim == 1 ) //single image // !decomposeStack must be first
        {
            id = addObject();
            MetaData::setValue(MDL_IMAGE, filename, id);
            MetaData::setValue(MDL_ENABLED, 1, id);
        }
        else //stack
        {
            FileName fnTemp;
            for (size_t i = 1; i <= image().ndim; ++i)
            {
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
    {
        REPORT_ERROR(ERR_IO_NOTEXIST, formatString("MetaDataDb::read: File doesn't exists: %s", inFile.c_str()) );
    }

    bool useCommentAsImage = false;
    this->_inFile = inFile;
    bool oldFormat=true;

    is.seekg(0, std::ios::beg);//reset the stream position to the beginning to start parsing

    if (line.find(FileNameVersion) != String::npos ||
        extFile == "xmd" || extFile == "star")
    {
        oldFormat = false;
        _comment.clear();

        // Skip comment parsing if we found the data key in the first line
        if (line.find("data_") != 0)
        {
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
            setComment(_comment);
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

        while (nextBlock(buffer, block))
            //startingPoint, remainingSize, firstData, secondData, firstloop))
        {
            BLOCK_NAME(block, blockName);
            if (blockRegExp.size() == 0 || regexec(&re, blockName.c_str(), (size_t) 0, NULL, 0)==0)
            {
                //Read column labels from the datablock that starts at firstData
                //Label ends at firstloop
                if ((_isColumnFormat = (block.loop != NULL)))
                {
                    _readColumnsStar(block, columnValues, desiredLabels, firstBlock);
                    // If block is empty, makes block.loop and block.end equal
                    if(block.loop == (block.end + 1))
                        block.loop--;
                    _readRowsStar(block, columnValues, desiredLabels);
                }
                else
                {
                    id = addObject();
                    _parsedLines = 1;
                    _readColumnsStar(block, columnValues, desiredLabels, firstBlock, id);
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
    }
    else if (line.find("Headerinfo columns:") != String::npos)
    {
        //This looks like an old DocFile, parse header
        std::cerr << "WARNING: ** You are using an old file format (DOCFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        is.ignore(256, ':'); //ignore all until ':' to start parsing column labels
        getline(is, line);
        ss.str(line);
        columnValues.push_back(new MDObject(MDL_UNDEFINED));
        columnValues.push_back(new MDObject(MDL_UNDEFINED));

        addLabel(MDL_IMAGE);
        _readColumns(ss, columnValues, desiredLabels);
        useCommentAsImage = true;
    }
    else
    {
        std::cerr << "WARNING: ** You are using an old file format (SELFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        //I will assume that is an old SelFile, so only need to add two columns
        columnValues.push_back(new MDObject(MDL_IMAGE));//addLabel(MDL_IMAGE);
        columnValues.push_back(new MDObject(MDL_ENABLED));//addLabel(MDL_ENABLED);
    }

    if (oldFormat)
        _readRows(is, columnValues, useCommentAsImage);

    //free memory of column values
    int nCols = columnValues.size();
    for (int i = 0; i < nCols; ++i)
        delete columnValues[i];

    is.close();
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
    {
        const char * srcName = MDL::label2Str(labelSrc).c_str();
        REPORT_ERROR(ERR_ARG_MISSING, formatString("Source label: '%s' doesn't exist on metadata", srcName));
    }
    md.addLabel(labelDest);
    std::vector<MDObject> values;
    getColumnValues(labelSrc, values);
    md.setColumnValues(values);
}

void MetaDataDb::renameColumn(MDLabel oldLabel, MDLabel newLabel)
{
    if (!containsLabel(oldLabel))
    {
        const char * srcName = MDL::label2Str(oldLabel).c_str();
        REPORT_ERROR(ERR_ARG_MISSING, formatString("Source label: '%s' doesn't exist on metadata",
                     srcName));
    }
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
    double max;
    MDObject result(column);
    aggregateSingle(result, AGGR_MAX, column);
    result.getValue(max);
    return max;

}

double MetaDataDb::getColumnMin(MDLabel column)
{
    double min;
    MDObject result(column);
    aggregateSingle(result, AGGR_MIN, column);
    result.getValue(min);
    return min;

}


void MetaDataDb::aggregateSingleInt(MDObject &mdValueOut, AggregateOperation op,
                                  MDLabel aggregateLabel)

{
    size_t aux = myMDSql->aggregateSingleSizeT(op,aggregateLabel);
    int aux2 = (int) aux;
    mdValueOut.setValue(aux2);
}


bool MetaDataDb::nextBlock(mdBuffer &buffer, mdBlock &block)
{
    BLOCK_INIT(block);
    if (buffer.size == 0)
        return false;
    // Search for data_ after a newline
    block.begin = BUFFER_FIND(buffer, "data_", 5);

    if (block.begin) // data_ FOUND!!!
    {
        block.begin += 5; //Shift data_
        size_t n = block.begin - buffer.begin;
        BUFFER_MOVE(buffer, n);
        //Search for the end of line
        char *newLine = BUFFER_FIND(buffer, "\n", 1);
        //Calculate length of block name, counting after data_
        block.nameSize = newLine - buffer.begin;
        //Search for next block if exists one
        //use assign and check if not NULL at same time
        if (!(block.end = BUFFER_FIND(buffer, "\ndata_", 6))) {
            block.end = block.begin + buffer.size; }
        else {
            block.end += 1; // to include terminal \n
        }
        block.loop = BUFFER_FIND(buffer, "\nloop_", 6);
        //If loop_ is not found or is found outside block
        //scope, the block is in column format
        if (block.loop)
        {
            if (block.loop < block.end)
                block.loop += 6; // Shift \nloop_
            else
                block.loop = NULL;
        }
        //Move buffer to end of block
        n = block.end - buffer.begin;
        BUFFER_MOVE(buffer, n);
        return true;
    }

    return false;
}

void MetaDataDb::aggregate(const MetaDataDb &mdIn, AggregateOperation op,
                         MDLabel aggregateLabel, MDLabel operateLabel, MDLabel resultLabel)
{
    std::vector<MDLabel> labels(2);
    std::vector<MDLabel> operateLabels(1);
    labels[0] = aggregateLabel;
    labels[1] = resultLabel;
    operateLabels[0]=operateLabel;
    init(&labels);
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
    init(&resultLabels);
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
    labels.push_back(resultLabel);
    init(&labels);
    mdIn.myMDSql->aggregateMdGroupBy(this, op, groupByLabels, operateLabel, resultLabel);
}

//-------------Set Operations ----------------------
void MetaDataDb::_setOperates(const MetaDataDb &mdIn,
                            const MDLabel label,
                            SetOperation operation)
{
    std::vector<MDLabel> labels;
    labels.push_back(label);
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
    labels.push_back(label);
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

void MetaDataDb::removeDisabled()
{
    if (containsLabel(MDL_ENABLED))
        removeObjects(MDValueLE(MDL_ENABLED, 0)); // Remove values -1 and 0 on MDL_ENABLED label
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
    labelsLeft.push_back(labelLeft);
    labelsRight.push_back(labelRight);
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
        init(&(MDin._activeLabels));
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
        init(&(MDin._activeLabels));
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
    init(&(mdIn._activeLabels));
    copyInfo(mdIn);
    mdIn.myMDSql->copyObjects(this, new MDQuery(n_images, first, sortLabel));
}

void MetaDataDb::selectSplitPart(const MetaDataDb &mdIn, size_t n, size_t part, const MDLabel sortLabel)
{
    size_t mdSize = mdIn.size();
    if (n > mdSize)
        REPORT_ERROR(ERR_MD, "selectSplitPart: Couldn't split a metadata in more parts than its size");
    if (part < 0 || part >= n)
        REPORT_ERROR(ERR_MD, "selectSplitPart: 'part' should be between 0 and n-1");
    _selectSplitPart(mdIn, n, part, mdSize, sortLabel);

}

void MetaDataDb::selectRandomSubset(const MetaDataDb &mdIn, size_t numberOfObjects, const MDLabel sortLabel)
{
    clear();

    MetaDataDb mdAux, mdAux2;
    mdAux.randomize(mdIn);
    mdAux2.selectPart(mdAux, 0, numberOfObjects);
    sort(mdAux2,sortLabel);
}

void MetaDataDb::selectPart(const MetaDataDb &mdIn, size_t startPosition, size_t numberOfObjects,
                          const MDLabel sortLabel)
{
    size_t mdSize = mdIn.size();
    if (startPosition < 0 || startPosition >= mdSize)
        REPORT_ERROR(ERR_MD, "selectPart: 'startPosition' should be between 0 and size()-1");
    init(&(mdIn._activeLabels));
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
        vd.push_back(dynamic_cast<const MDRowSql&>(row));
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
