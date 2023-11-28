/***************************************************************************
 *
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Jan Horacek (xhorace4@fi.muni.cz)
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

/**
 * This file defines abstract class MetaData which defines an API for all
 * MetaData implementations.
 */

#ifndef CORE_METADATA_H
#define CORE_METADATA_H

#include <cstddef>
#include <map>
#include <cmath>
#include <memory>
#include <iostream>

#include "xmipp_error.h"
#include "xmipp_filename.h"
#include "metadata_label.h"
#include "metadata_writemode.h"
#include "metadata_base_it.h"
#include "metadata_static.h"
#include "metadata_query.h"
#include "choose.h"

/** @defgroup MetaData Metadata Stuff
 * @ingroup DataLibrary
 * @{
 */

#define BAD_OBJID 0
#define BAD_INDEX -1

#define FILENAME_XMIPP_STAR "# XMIPP_STAR_1"
#define FILENAME_XMIPP_SQLITE "SQLite format 3"
#define DEFAULT_BLOCK_NAME "noname"

// FIXME: deprecated
// Preffered iterating is on right side of these macros
#define FOR_ALL_OBJECTS_IN_METADATA(__md) for (size_t objId : __md.ids())

// FIXME: deprecated
#define FOR_ALL_ROWS_IN_METADATA(__md) for (auto& row : __md)

#define END_OF_LINE() ((char*) memchr (iter, '\n', end-iter))

/** Which are the blocks available in a metadata */
void getBlocksInMetaDataFile(const FileName &inFile, StringVector& blockList);
bool existsBlockInMetaDataFile(const FileName &inFile, const String& inBlock);
bool existsBlockInMetaDataFile(const FileName &inFileWithBlock);

class MDValueGenerator;

/** Struct to hold a char * pointer and a size
 * this will be useful for parsing metadata
 */
typedef struct {
    char * begin;
    size_t size;
}
mdBuffer;

/// Some macros to use the buffer
#define BUFFER_CREATE(b) mdBuffer b; b.begin = nullptr; b.size = 0
#define BUFFER_COPY(b1, b2) mdBuffer b2; b2.begin = b1.begin; b2.size = b1.size
#define BUFFER_MOVE(b, n) b.begin += n; b.size -= n
#define BUFFER_FIND(b, str, n) (char*) _memmem(b.begin, b.size, str, n)

typedef struct {
    char * begin; //Position of _dataXXX on buffer
    size_t nameSize; //Number of charater of block name, counting after _data
    char * end; //Position just before next _dataXXX or end of buffer
    char * loop; //Position of _loop if exists, NULL otherwise
}
mdBlock;
/// Some macros to use the block pointers
#define BLOCK_CREATE(b) mdBlock b; b.begin = b.end = b.loop = nullptr; b.nameSize = 0
#define BLOCK_INIT(b) b.begin = b.end = b.loop = nullptr; b.nameSize = 0
#define BLOCK_NAME(b, s) s.assign(b.begin, b.nameSize)

class ObjectDoesNotExist: public std::logic_error {
public:
    ObjectDoesNotExist(MDLabel label, size_t id) : std::logic_error("Label '" + MDL::label2Str(label) + " does not exist for id " + std::to_string(id)) {};
    ObjectDoesNotExist(size_t id, const FileName &fn) : std::logic_error("Row with id " + std::to_string(id) + " does not exist" + (fn.isEmpty() ? "" : " in " + fn)) {};
};

class NotImplemented : public std::logic_error {
public:
    NotImplemented(const std::string& str) : std::logic_error(str) {};
};

class ColumnDoesNotExist : public std::logic_error {
public:
    ColumnDoesNotExist(MDLabel label, const FileName &fn) : std::logic_error("Column '" + MDL::label2Str(label) + " does not exist" + (fn.isEmpty() ? "" : " in " + fn)) {};
    ColumnDoesNotExist(const std::string &msg) : std::logic_error(msg) {};
};


/**
 * Definition of API of all MetaDatas.
 *
 * ### General information about metadata
 *
 * MetaData* class(es) implement a way to store metadata in xmipp programs.
 * - MetaData is represented as as a database-like table.
 * - Each intersection of row and column contains MDObject instance.
 * - Each row has its label MDLabel.
 * - Columns of all rows are same.
 * - Metadata could be loaded from file, stored to file, iterated over, rows
 *   can be added, removed, changed etc.
 *
 * There was a single original database (MetaDataDb) implementation of MetaData
 * till 2021, when MetaData were split into MetaDataDb & MetaDataVec
 * with aim to achive higner speeds via saving metadata in std::vectors instead
 * of sql database.
 *
 * Current MetaData implementation:
 *  1. metadata_base.(h|cpp): common MetaData API definition
 *  2. metadata_vec.(h|cpp): vector MetaData implementation
 *  3. metadata_db.(h|cpp): old datababse MetaData implementation
 *
 * ### MetaDataBase
 *
 * MetaData class cannot be instantiated however it could be passed to functions
 * as reference or pointer where only common MetaData functions are required:
 *  - `void foo(MetaData&);`
 *  - `void boo(MetaData*);`
 *
 * This code doesn't compile:
 *  `void roo(MetaData& md) { md.someDbSpeficicOperation(); }`
 *
 * This code compiles:
 *  `void roo(MetaDataDb& md) { md.someDbSpeficicOperation(); }`
 *
 * ### Iterating over MetaData
 *
 * You can iterate directly over abstract MetaData as well as over specific MetaData*.
 * You can iterate over:
 *  1. Ids of rows (arbitrary size_t): `for (size_t id : md.ids())`
 *  2. Rows: `for (const MDRow& row : md)`
 *
 * ### Useful information
 *
 * - See metadata.h file for general info about Metadata.
 * - Never rely on any properties of IDs! That is, do NOT assume that IDs are
 *   e.g. continuous or that they are ascending.
 * - There is a hierarchy of MDRow classes copying MetaData hierarchy: MDRow, MDRowVec, MDRowSql.
 *   If you need row-specific function, do a dynamic cast:
 *   `for (const MDRow& row : md) { const MDRowVec& rowv = dynamic_cast<MDRowVec&>(row); }`
 *   You need to make sure type of md is MetaDataVec (otherwise exception is thrown).
 */
class MetaData {
private:
    bool _isColumnFormat; // Format for the file, column or row formatted
protected:
    /* Allows a fast search for pairs where the value is
     * a string, i.e. looking for filenames which is quite usual.
     */
    std::map<String, size_t> _fastStringSearch;
    MDLabel _fastStringSearchLabel;
    String _comment; ///< A general comment for the MetaData file
    ///comment is wraped in char_max length lines
#define line_max 70

    
    int _precision = 1000;

    /**Input file name
     * Where does this MetaData come from/go to be stored?
     */
    FileName _inFile;

    /** This two variables will be used to read the metadata information (labels and size)
     * or maybe a few rows only
     */
    size_t _maxRows = 0, _parsedLines = 0;

    void copyInfo(const MetaData& md);
    double precision() const;

    virtual void readStar(const FileName &filename, const std::vector<MDLabel> *desiredLabels,
                          const String &blockRegExp, bool decomposeStack);

    virtual void _readColumns(std::istream& is, std::vector<MDObject*> & columnValues,
                              const std::vector<MDLabel>* desiredLabels = nullptr);

    virtual void _readColumnsStar(mdBlock &block,
                                  std::vector<MDObject*> & columnValues,
                                  const std::vector<MDLabel>* desiredLabels,
                                  bool addColumns = true,
                                  size_t id = BAD_OBJID);

    /* Helper function to parse an MDObject and set its value.
     * The parsing will be from an input stream(istream)
     * and if parsing fails, an error will be raised
     */
    virtual void _parseObject(std::istream &is, MDObject &object, size_t id = BAD_OBJID);

    virtual void _parseObjects(std::istream &is, std::vector<MDObject*> & columnValues,
                               const std::vector<MDLabel> *desiredLabels, bool firstTime) = 0;

    /* This function will be used to parse the rows data
     * having read the columns labels before and setting which are desired
     * the useCommentAsImage is for compatibility with old DocFile format
     * where the image were in comments
     */
    virtual void _readRows(std::istream& is, std::vector<MDObject*>& columnValues, bool useCommentAsImage);

    /** This function will be used to parse the rows data in START format
     * @param[out] columnValues MDRow with values to fill in
     * @param pchStart pointer to the position of '_loop' in memory
     * @param pEnd  pointer to the position of the next '_data' in memory
     * @param maxRows if this number if greater than 0, only this number of rows will be parsed.
     */
    virtual void _readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues,
                               const std::vector<MDLabel> *desiredLabels);

public:
    /** Filename used in the read command, useful to write Error messages
     *
     */
    bool isMetadataFile;

    /** @name Constructors
     *  @{
     */

    /** Empty Constructor.
     *
     * The MetaData is created with no data stored on it. You can fill in it programmatically
     * or by a later reading from a MetaData file or old Xmipp formatted type.
     * if labels vectors is passed this labels are created on metadata
     */
    MetaData() = default;
    virtual ~MetaData();

    /**Clear all data
     */
    virtual void clear();
    /** @} */

    /** @name Getters and setters
     * @{
     */

    /**Return true if the metadata is in column format.
     */
    virtual bool isColumnFormat() const { return _isColumnFormat; }

    /** Prevent from parsing all rows from the metadata.
     * When reading from file, only maxRows will be read.
     */
    virtual void setMaxRows(size_t maxRows=0) { _maxRows = maxRows; }

    /** Return the number of lines in the metadata file.
     * Serves to know the number of items even is read with
     * maxRows != 0
     */
    virtual size_t getParsedLines() { return _parsedLines; }

    /**Set precision (number of decimal digits) use by operator == when comparing
     * metadatas with double data. "2" is a good value for angles
     */
    virtual void setPrecission(int _precision) { this->_precision = (int)pow (10,_precision); }

    /** Set to false for row format (parameter files).
     *  set to true  for column format (this is the default) (docfiles)
     */
    virtual void setColumnFormat(bool column) { _isColumnFormat = column; }

    /** Export medatada to xml file.
     *
     */
    virtual void writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const = 0;

    /** Write metadata in text file as plain data without header.
     *
     */
    virtual void writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const = 0;

    virtual void writeStar(const FileName &outFile, const String & blockName, WriteModeMetaData mode) const;

    /**Get Header Comment.
     * the comment will appear in second line.
     */
    virtual String getComment() const { return this->_comment; }

    /**Set Header Comment.
     * the comment will appear in second line
     */
    virtual void setComment(const String &newComment = "No comment") { this->_comment = newComment; }

    /**Get metadata filename.
     */
    virtual FileName getFilename() const { return this->_inFile; }

    /**Set metadata filename.
     */
    virtual void setFilename(const FileName &_filename) { this->_inFile = _filename; }

    /**Get safe access to active labels.
     */
    virtual std::vector<MDLabel> getActiveLabels() const = 0;

    /**Get maximum string length of column values.
    */
    virtual int getMaxStringLength(const MDLabel thisLabel) const = 0;

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */

    /** Set the value of all objects in an specified column (both value and column are specified in mdValueIn)
    */
    virtual bool setValueCol(const MDObject &mdValueIn) = 0;

    /**Set the value of all objects in an specified column.
     * @code
     * MetaData md;
     * md.setValueCol(MDL_IMAGE, "images/image00011.xmp");
     * @endcode
     */
    template<class T>
    bool setValueCol(const MDLabel label, const T &valueIn) {
        return setValueCol(MDObject(label, valueIn));
    }

    /** Set the value for some label.
     * to the object that has id 'objectId'
     * or to 'activeObject' if is objectId=-1.
     * This is one of the most used functions to programatically
     * fill a metadata.
     * @code
     * MetaData md;
     * size_t id = md.addObject();
     * md.setValue(MDL_IMAGE, "images/image00011.xmp",id);
     * md.setValue(MDL_ANGLE_ROT, 0.,id);
     * @endcode
     */
    template<class T>
    bool setValue(const MDLabel label, const T &valueIn, size_t id) {
        return setValue(MDObject(label, valueIn), id);
    }

    virtual size_t addRow(const MDRow &row) = 0;

    MetaData& operator=(const MetaData &md) {
        if (&md == this) {
            return *this;
        }
        this->clear();
        this->copyInfo(md);
        for (const auto& row : md)
            this->addRow(row);
        return *this;
    };

    /** This functions are using MDObject for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */
    virtual bool setValue(const MDObject &mdValueIn, size_t id) = 0;
    virtual bool getValue(MDObject &mdValueOut, size_t id) const = 0;
    // FIXME: deprecated; this function should return MDObject, it is const-unsafe

    /** Get the value of some label.
     * from the object that has id 'objectId'
     * or from 'activeObject' if objectId=-1.
     * @code
     * MetaData md;
     * md.read("images.xmd");
     * FileName imageFn;     *
     * FOR_ALL_OBJECTS_IN_METADATA(md)
     * {
     *      md.getValue(MDL_IMAGE, imageFn);
     *      std::out << "Image: " << imageFn);
     * }
     * @endcode
     */
    template<class T>
    T getValue(const MDLabel label, size_t id) const {
        MDObject mdValueOut(label);
        if (!getValue(mdValueOut, id))
            throw ObjectDoesNotExist(label, id);
        return mdValueOut.getValue2(T());
    }

    template<class T>
    bool getValue(const MDLabel label, T &valueOut, size_t id) const {
        MDObject mdValueOut(label);
        if (!getValue(mdValueOut, id))
            return false;
        valueOut = mdValueOut.getValue2(T());
        return true;
    }

    template<class T>
    const T& getValueOrAbort(const MDLabel label, size_t id) const {
        return getValue<T>(label, id);
    }

    template<class T>
    T& getValueOrAbort(const MDLabel label, size_t id) {
        return getValue<T>(label, id);
    }

    template<class T>
    void getValueOrAbort(const MDLabel label, T &valueOut, size_t id) const { // FIXME: deprecated
        if (!getValue(label, valueOut,id))
            REPORT_ERROR(ERR_ARG_MISSING,(String)"Cannot find label: " + MDL::label2Str(label));
    }

    template <typename T>
    const T getValueOrDefault(const MDLabel label, size_t id, const T &_default) const {
        MDObject mdValueOut(label);
        if (!getValue(mdValueOut, id))
            return _default;
        return mdValueOut.getValue2(T());
    }

    template <typename T>
    T getValueOrDefault(const MDLabel label, size_t id, T &_default) {
        MDObject mdValueOut(label);
        if (!getValue(mdValueOut, id))
            return _default;
        return mdValueOut.getValue2(T());
    }

    template <typename T, typename T1>
    void getValueOrDefault(const MDLabel label, T &valueOut, size_t id, const T1 &_default) const { // FIXME: deprecated
        if (!getValue(label, valueOut,id))
            valueOut = (T) _default;
    }

    /** Get all values of a column as a vector.
     */
    template<class T>
    std::vector<T> getColumnValues(const MDLabel label) const {
        std::vector<T> result;
        MDObject mdValueOut(label);
        std::vector<size_t> objectsId;
        findObjects(objectsId);
        size_t n = objectsId.size();
        result.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            if (!getValue(mdValueOut, objectsId[i]))
                throw std::logic_error("Column does not exist!");
            result.emplace_back(mdValueOut.getValue2(T()));
        }
        return result;
    }

    // FIXME: deprecated
    template<class T>
    void getColumnValues(const MDLabel label, std::vector<T>& out) const {
        out = this->getColumnValues<T>(label);
    }

    // FIXME: add proper return value
    virtual bool getRowValues(size_t id, std::vector<MDObject> &values) const = 0;

    /** Get all values of a column as a vector.
     */
    // FIXME: add proper return value
    virtual void getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const = 0;

    /** Set all values of a column as a vector.
     * The input vector must have the same size as the Metadata.
     */
    template<class T>
    void setColumnValues(const MDLabel label, const std::vector<T> &valuesIn) {
        if (valuesIn.size() != size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Input vector must be of the same size as the metadata");
        size_t n = 0;
        for (size_t objId : this->ids())
            this->setValue(label, valuesIn[n++], objId);
    }

    virtual void setColumnValues(const std::vector<MDObject> &valuesIn) = 0;

    virtual std::unique_ptr<MDRow> getRow(size_t id) = 0;
    virtual std::unique_ptr<const MDRow> getRow(size_t id) const = 0;

    /** Set label values from string representation.
     */
    virtual bool setValueFromStr(const MDLabel label, const String &value, size_t id);

    /** Get string representation from label value.
     */
    virtual bool getStrFromValue(const MDLabel label, String &strOut, size_t id) const;

    /**Check whether the metadata is empty.
     */
    virtual bool isEmpty() const { return size() == 0; }

    /**Number of objects contained in the metadata.
     */
    virtual size_t size() const = 0;

    /** Check whether a label is contained in metadata.
     */
    virtual bool containsLabel(const MDLabel label) const = 0;

    /** Add a new label to the metadata.
     * By default the label is added at the end,
     * if the position is specified and is between 0 and n-1
     * the new label is inserted at that position.
     */
    virtual bool addLabel(const MDLabel label, int pos = -1) = 0;

    /** Remove a label from the metadata.
     * The data is still in the table. If you want to remove the data,
     * make a copy of the MetaData.
     */
    virtual bool removeLabel(const MDLabel label) = 0;

    /** Remove all the labels from the metadata but the
     * ones given in labels vector.
     */
    void keepLabels(const std::vector<MDLabel> &labels);

    /** Adds a new, empty object to the objects map. If objectId == -1
     * the new ID will be that for the last object inserted + 1, else
     * the given objectId is used. If there is already an object whose
     * objectId == input objectId, just removes it and creates an empty
     * one
     */
    virtual size_t addObject() = 0;

    /** Import objects from another metadata.
     * @code
     * //Import object 1000 from metadata B into metadata A
     * A.importObject(B, 1000);
     * //Import all objects with rotational angle greater that 60
     * A.importObjects(B, MDValuesGT(MDL_ANGLE_ROT, 60));
     * //Import all objects
     * A.importObjects(B);
     * @endcode
     */
    virtual void importObject(const MetaData &md, const size_t id, bool doClear=true) = 0;
    virtual void importObjects(const MetaData &md, const std::vector<size_t> &objectsToAdd, bool doClear=true) = 0;
    virtual void importObjects(const MetaData &md, const MDQuery &query, bool doClear=true) = 0;

    /** Remove the object with this id.
    * Returns true if the object was removed or false if
    * the object did not exist
    */
    virtual bool removeObject(size_t id) = 0;

    /** Removes the collection of objects of given vector id's
     * NOTE: The iterator will point to the first object after any of these
     * operations
     */
    virtual void removeObjects(const std::vector<size_t> &toRemove) = 0;

    /** Removes objects from metadata.
     * return the number of deleted rows
     * if not query, all objectes are removed
     * Queries can be used in the same way
     * as in the importObjects function
     */
    virtual int removeObjects(const MDQuery&) = 0;
    virtual int removeObjects() = 0;

    /** @} */

    /** @name Iteration functions
     * @{
     */

    /** Return the object id of the first element in metadata. */
    virtual size_t firstRowId() const = 0;
    virtual size_t firstObject(const MDQuery&) const = 0;

    /** Goto last metadata object.*/
    virtual size_t lastRowId() const = 0;

    /** @name Search operations
     * @{
     */

    /** Find all objects that match a query.
     * if called without query, all objects are returned
     * if limit is provided only return a maximun of 'limit'
     */
    virtual void findObjects(std::vector<size_t> &objectsOut, const MDQuery &query) const = 0;
    virtual void findObjects(std::vector<size_t> &objectsOut, int limit = -1) const = 0;

    virtual size_t countObjects(const MDQuery&) const = 0;
    virtual bool containsObject(size_t objectId) const = 0;
    virtual bool containsObject(const MDQuery&) const = 0;

    /** @} */

    /** @name I/O functions
     * @{
     */

    virtual void _writeRows(std::ostream &os) const = 0;

    virtual void write(const FileName &outFile, WriteModeMetaData mode=MD_OVERWRITE) const = 0;
    virtual void write(std::ostream &os, const String & blockName="",WriteModeMetaData mode=MD_OVERWRITE) const = 0;
    virtual void print() const { this->write(std::cout); }

    /** Append data lines to file.
     * This function can be used to add new data to
     * an existing metadata. Now should be used with
     * files with only one metadata, maybe can be extended later.
     * For now it will not check any compatibility beetween the
     * existent metadata and the new data to append.
     */
    void append(const FileName &outFile) const;

    /** Read data from file. Guess the blockname from the filename
     * @code
     * inFilename="first@md1.doc" -> filename = md1.doc, blockname = first
     * @endcode
     */
    virtual void read(const FileName &inFile, const std::vector<MDLabel> *desiredLabels = nullptr, bool decomposeStack=true) = 0;
    /** @} */

    /** @name Set Operations
     * @{
     */

    /**Remove rows with MDL_ENABLED = -1 if this label is present
     */
    virtual void removeDisabled();

    /** Select random subset */
    virtual void selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects,
                                    const MDLabel sortLabel=MDL_OBJID) = 0;

    /** Select some part from Metadata.
     * Select elements from input Metadata
     * at some starting position
     * if the numberOfObjects is -1, all objects
     * will be returned from startPosition to the end.
    */
    virtual void selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                            const MDLabel sortLabel=MDL_OBJID) = 0;

    /** @} */

    /** @name Iterators
     *  @{
     */
    friend struct MDBaseRowIterator<false>;
    friend struct MDBaseRowIterator<true>;

    /* To allow iteration over abstract MetaData, there is a hierarchy of iterators.
       - MDBaseRowIterator, MDVecRowIterator, MDDbRowIterator
       - MDBaseIdIterator, MDVecIdIterator, MDDbIdIterator
     * Plus there are rowIterator & idIterator, which are common for all instances
     * and just hold unique_ptr to appropriate implementation.
     */
    template <bool IsConst>
    struct rowIterator {
    private:
        std::unique_ptr<MDBaseRowIterator<IsConst>> impl;
    public:
        rowIterator(std::unique_ptr<MDBaseRowIterator<IsConst>> impl) : impl(std::move(impl)) {}
        rowIterator(rowIterator const& right) : impl(std::move(right.impl->clone())) {}
        rowIterator& operator=(rowIterator const& right) {
            impl = std::move(right.impl->clone());
            return *this;
        }
        rowIterator& operator++() {
            impl->increment();
            return *this;
        }
        bool operator==(const rowIterator<IsConst>& other) const { return *(other.impl) == *(this->impl); }
        bool operator!=(const rowIterator<IsConst>& other) const { return !(*this == other); }
        typename TypeHelpers::choose<IsConst, const MDRow&, MDRow&>::type operator*() const { return **impl; }
    };

    using iterator = rowIterator<false>;
    using const_iterator = rowIterator<true>;

    virtual iterator begin() = 0;
    virtual iterator end() = 0;

    virtual const_iterator begin() const = 0;
    virtual const_iterator end() const = 0;


    template <bool IsConst>
    struct idIterator {
    private:
        std::unique_ptr<MDBaseIdIterator<IsConst>> impl;
    public:
        idIterator(std::unique_ptr<MDBaseIdIterator<IsConst>> impl) : impl(std::move(impl)) {}
        idIterator(idIterator const& right) : impl(std::move(right.impl->clone())) {}
        idIterator& operator=(idIterator const& right) {
            impl = std::move(right.impl->clone());
            return *this;
        }
        idIterator& operator++() {
            impl->increment();
            return *this;
        }
        bool operator==(const idIterator<IsConst>& other) const { return *(other.impl) == *(this->impl); }
        bool operator!=(const idIterator<IsConst>& other) const { return !(*this == other); }
        size_t operator*() const { return **impl; }
    };

    using id_iterator = idIterator<false>;
    using id_const_iterator = idIterator<true>;

    // This proxy allows to implement non-const & const iterator in one templated class.
    template <bool IsConst>
    struct IdIteratorProxy {
        typename TypeHelpers::choose<IsConst, const MetaData&, MetaData&>::type _md;

        IdIteratorProxy(typename TypeHelpers::choose<IsConst, const MetaData&, MetaData&>::type md) : _md(md) { }
        typename TypeHelpers::choose<IsConst, MetaData::id_const_iterator, MetaData::id_iterator>::type begin() { return _md.id_begin(); };
        typename TypeHelpers::choose<IsConst, MetaData::id_const_iterator, MetaData::id_iterator>::type end() { return _md.id_end(); };
    };

    virtual id_iterator id_begin() = 0;
    virtual id_iterator id_end() = 0;

    virtual id_const_iterator id_begin() const = 0;
    virtual id_const_iterator id_end() const = 0;

    virtual IdIteratorProxy<false> ids() { return IdIteratorProxy<false>(*this); };
    virtual IdIteratorProxy<true> ids() const { return IdIteratorProxy<true>(*this); };

    /** @} */

    /** Fill column with constant value
     */
    virtual void fillConstant(MDLabel label, const String &value) = 0;

    /** Fill column with random value
     * mode should be: uniform, gaussian or student
     * op1, op2 and op2 are interpreted for each mode:
     * uniform: op1 and op2 are the limits of the interval
     * gaussian: op1 and op2 are mean and std
     * student: same as gaussian and use op3
     */
    virtual void fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3=0.) = 0;

    /** Fill lineal, starting at some value and with some step */
    virtual void fillLinear(MDLabel label, double initial, double step) = 0;

    /** Copy all values from one column to another.
     * Source column should exist
     */
    virtual void copyColumn(MDLabel labelDest, MDLabel labelSrc) = 0;

    /** Same as previous, but copy to another metadata */
    virtual void copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc) = 0;

    /** Rename column.
     *
     */
    virtual void renameColumn(MDLabel oldLabel, MDLabel newLabel) = 0;

    /** Rename several columns. This is an expensive operations so if several
     * columns need to be changed do it using this function instead one by one
     */
    virtual void renameColumn(const std::vector<MDLabel> &oldLabel,
            const std::vector<MDLabel> &newLabel) = 0;

    virtual bool nextBlock(mdBuffer &buffer, mdBlock &block);
};//class MetaData

/** print metadata
 *
 */
std::ostream& operator<<(std::ostream& o, const MetaData & mD);

/** @} */

/** Convert string to write mode metadata enum.
 *
 */
WriteModeMetaData metadataModeConvert (String mode);

bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label);

#endif
