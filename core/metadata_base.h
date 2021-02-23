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

#ifndef CORE_METADATA_H
#define CORE_METADATA_H

#include <cstddef>
#include <map>
#include <cmath>
#include <memory>

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


/** Iterate over all elements in MetaData
 *
 * This macro is used to generate loops over all elements in the MetaData.
 * At each iteration the 'active object' is changed so you can perform
 * the set and get default method on the MetaData.
 *
 * @code
 * MetaData md;
 *   //...
 * FOR_ALL_OBJECTS_IN_METADATA(md)
 * {
 *     String imageFile;
 *     md.getValue(MDL_IMAGE, imageFile);
 *     std::cout << "Image file: " << imageFile << " ";
 * }
 * @endcode
 */

#define FOR_ALL_OBJECTS_IN_METADATA(__md)   for (auto& row : __md)
// auto due to const/non-const

/* #define FOR_ALL_ROWS_IN_METADATA(__md) \
 *         for(MDRowIterator __iter(__md); __iter.hasNext(); __iter.moveNext(__md))
 */

/** Iterate over all elements of two MetaData at same time.
 *
 * This macro is useful to iterate over two MetaData with the same
 * number of elements and performs operations to elements in both of them.
 * At each iteration the 'active objects' in both MetaData are changed.
 *
 * @code
 * MetaData mdA, mdB, mdC;
 *  //Iterate over MetaData mdA and mdB
 *  //take image from the first and tilt angle from the second
 *  //and create a new MetaData.
 * FOR_ALL_OBJECTS_IN_METADATA2(mdA, mdB)
 * {
 *     String imageFile;
 *     double angle;
 *     mdA.getValue(MDL_IMAGE, imageFile,__iter.objId);
 *     mdB.getValue(MDL_ANGLE_TILT, angle,__iter2.objId);
 *     size_t objId=mdC.addObject();
 *     mdC.setValue(MDL_IMAGE, imageFile,objId);
 *     mdC.setValue(MDL_ANGLE_TILT, angle,objId);
 * }
 * @endcode
 */

#define FOR_ALL_OBJECTS_IN_METADATA2(__md, __md2) \
        for(MDIterator __iter(__md), __iter2(__md2);\
             __iter.hasNext() && __iter2.hasNext(); \
             __iter.moveNext(), __iter2.moveNext())

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


/** Class to manage data files.
 *
 * The MetaData class manages all procedures related to
 * metadata. MetaData is intended to group together old
 * Xmipp specific files like Docfiles, Selfiles, etc..
 *
 */
class MetaData {
public:
    // Allows a fast search for pairs where the value is
    // a string, i.e. looking for filenames which is quite
    // usual
    std::map<String, size_t> _fastStringSearch;
    MDLabel _fastStringSearchLabel;
    String _path; ///< A parameter stored on MetaData Files
    String _comment; ///< A general comment for the MetaData file
    ///comment is wraped in char_max length lines
#define line_max 70

    bool _isColumnFormat; ///< Format for the file, column or row formatted
    int _precision;
    /**Input file name
     * Where does this MetaData come from/go to be stored?
     */
    FileName _inFile;

    /** What labels have been read from a docfile/metadata file
     * and/or will be stored on a new metadata file when "save" is
     * called
     **/
    std::vector<MDLabel> _activeLabels;

    /** When reading a column formatted file, if a label is found that
     * does not exist as a MDLabel, it is ignored. For further
     * file processing, such columns must be ignored and this structure
     * allows to do that
     **/
    std::vector<unsigned int> _ignoreLabels;

    /** This two variables will be used to read the metadata information (labels and size)
     * or maybe a few rows only
     */
    size_t _maxRows, _parsedLines;

    void copyInfo(const MetaData& md);

public:
    /** @name Constructors
     *  @{
     */

    /** Empty Constructor.
     *
     * The MetaData is created with no data stored on it. You can fill in it programmatically
     * or by a later reading from a MetaData file or old Xmipp formatted type.
     * if labels vectors is passed this labels are created on metadata
     */
    MetaData();
    MetaData(const std::vector<MDLabel> *labelsVector);

    /** From File Constructor.
     *
     * The MetaData is created and data is read from provided FileName. Optionally, a vector
     * of labels can be provided to read just those required labels
     */
    MetaData(const FileName &fileName, const std::vector<MDLabel> *desiredLabels = nullptr);

    /** Copy constructor
     *
     * Created a new metadata by copying all data from an existing MetaData object.
     */
    MetaData(const MetaData &md);

    /** Assignment operator
     *
     * Copies MetaData from an existing MetaData object.
     */
    virtual MetaData& operator=(const MetaData &md);

    /** Destructor
     *
     * Frees all used memory and destroys object.
     */
    ~MetaData();

    /**Clear all data
     */
    virtual void clear() = 0;
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

    bool nextBlock(mdBuffer &buffer, mdBlock &block); // TODO

    /** Export medatada to xml file.
     *
     */
    virtual void writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const = 0;

    /** Write metadata in text file as plain data without header.
     *
     */
    virtual void writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const = 0;

    /**Get path.
     */
    virtual String getPath() const { return this->_path; }

    /**Set Path.
     * the path will appear in first line
     */
    virtual void setPath(const String &newPath = "") {
        const size_t length = 512;
        char _buffer[length];
        this->_path = (newPath == "") ? String(getcwd(_buffer, length)) : newPath;
    }

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
    virtual std::vector<MDLabel> getActiveLabels() const { return this->_activeLabels; }

    /**Get unsafe pointer to active labels.
     */
    virtual std::vector<MDLabel> *getActiveLabelsAddress() const {
        return (std::vector<MDLabel>*) (&this->_activeLabels);
    }

    /**Get maximum string length of column values.
    */
    virtual int getMaxStringLength(const MDLabel thisLabel) const = 0;

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */


    /** Set the value of all objects in an specified column (both value and column are specified in mdValueIn)
    */
    virtual bool setValueCol(const MDObject &mdValueIn);

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

    /** This functions are using MDObject for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */
    virtual bool setValue(const MDObject &mdValueIn, size_t id);
    virtual bool getValue(MDObject &mdValueOut, size_t id) const;
    /** Filename used in the read command, useful to write Error messages
     *
     */
    FileName eFilename; // TODO
    String extFile; //Filename extension TODO
    bool isMetadataFile; // TODO

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
    bool getValue(const MDLabel label, T &valueOut, size_t id) const {
        MDObject mdValueOut(label);
        if (!getValue(mdValueOut, id))
            return false;
        mdValueOut.getValue(valueOut);
        return true;
    }

    template<class T>
    void getValueOrAbort(const MDLabel label, T &valueOut, size_t id) const {
        if (!getValue(label, valueOut,id))
            REPORT_ERROR(ERR_ARG_MISSING,(String)"Cannot find label: " + MDL::label2Str(label));
    }

    template <typename T, typename T1>
    void getValueOrDefault(const MDLabel label, T &valueOut, size_t id, const T1 &_default) const {
        if (!getValue(label, valueOut,id))
            valueOut = (T) _default;
    }

    /** Get all values of a column as a vector.
     */
    template<class T>
    void getColumnValues(const MDLabel label, std::vector<T> &valuesOut) const {
        T value;
        MDObject mdValueOut(label);
        std::vector<size_t> objectsId;
        findObjects(objectsId);
        size_t n = objectsId.size();
        valuesOut.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            getValue(mdValueOut, objectsId[i]);
            mdValueOut.getValue(value);
            valuesOut[i] = value;
        }
    }

    virtual bool getRowValues(size_t id, std::vector<MDObject> &values) const = 0;

    /** Get all values of a column as a vector.
     */
    virtual void getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const = 0;

    /** Set all values of a column as a vector.
     * The input vector must have the same size as the Metadata.
     */
    template<class T>
    void setColumnValues(const MDLabel label, const std::vector<T> &valuesIn) {
        bool addObjects=false;
        if (size()==0)
            addObjects=true;
        if (valuesIn.size()!=size() && !addObjects)
            REPORT_ERROR(ERR_MD_OBJECTNUMBER, "Input vector must be of the same size as the metadata");
        if (!addObjects) {
            size_t n = 0;
            for (auto row : *this)
                row.setValue(label, valuesIn[n++]);
        } else {
            size_t nmax=valuesIn.size();
            for (size_t n=0; n<nmax; ++n)
                setValue(label, valuesIn[n], addObject());
        }
    }

    virtual void setColumnValues(const std::vector<MDObject> &valuesIn) = 0;

    /*bool bindValue(size_t id) const;

    bool initGetRow(bool addWhereClause) const;
    bool execGetRow(MDRow &row) const;
    void finalizeGetRow(void) const;
    bool getRow(MDRow &row, size_t id) const;
    bool getAllRows(std::vector<MDRow> &rows) const;
    bool getRow2(MDRow &row, size_t id) const;

    ** Copy all the values in the input row in the current metadata
    bool initSetRow(const MDRow &row);
    bool execSetRow(const MDRow &row, size_t id);
    void finalizeSetRow(void);
    bool setRow(const MDRow &row, size_t id);
    bool setRow2(const MDRow &row, size_t id);

    ** Add a new Row and set values, return the objId of newly added object
    bool initAddRow(const MDRow &row);
    bool execAddRow(const MDRow &row);
    void finalizeAddRow(void);
    size_t addRow(const MDRow &row);
    void addRowOpt(const MDRow &row);
    void addRows(const std::vector<MDRow> &rows);
    void addMissingLabels(const MDRow &row);
    size_t addRow2(const MDRow &row);*/

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
    virtual bool containsLabel(const MDLabel label) const {
        return vectorContainsLabel(this->_activeLabels, label);
    }

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
    virtual bool keepLabels(const std::vector<MDLabel> &labels) = 0;

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
     * A.importObjects(B);     *
     * @endcode
     */
    virtual void importObject(const MetaData &md, const size_t id, bool doClear=true) = 0;
    virtual void importObjects(const MetaData &md, const std::vector<size_t> &objectsToAdd, bool doClear=true) = 0;

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

    /** Add item id.
     * From 1 to last.
     */
    // void addItemId();

    /** Remove item id.*/
    // void removeItemId();

    /** @} */

    /** @name Iteration functions
     * @{
     */

    /** Return the object id of the first element in metadata. */
    virtual size_t firstRowId() const;
    virtual size_t firstObject() const { return firstRowId(); } // deprecated: bad name
    virtual size_t firstObject(const MDQuery&) const = 0;

    /** Goto last metadata object.*/
    virtual size_t lastRowId() const;
    virtual size_t lastObject() const { return lastRowId(); } // deprecated: bad name

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

    virtual void write(const FileName &outFile, WriteModeMetaData mode=MD_OVERWRITE) const = 0;
    virtual void write(std::ostream &os, const String & blockName="",WriteModeMetaData mode=MD_OVERWRITE) const = 0;
    virtual void print() const = 0;

    /** Append data lines to file.
     * This function can be used to add new data to
     * an existing metadata. Now should be used with
     * files with only one metadata, maybe can be extended later.
     * For now it will not check any compatibility beetween the
     * existent metadata and the new data to append.
     */
    // void append(const FileName &outFile) const;

    /** Check if block exists in metadata file
     * input full parh block@filename
     * return false if metadata block does not exits
     */
    // bool existsBlock(const FileName &_inFile);

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

    /** Returns Max and Min values from a column in metadata
     * These functions can only be used for labels of type double
     */
    // double getColumnMax(MDLabel column);

    // double getColumnMin(MDLabel column);

    /** Basic operations on columns data.
     * Mainly perform replacements on string values and
     * basic algebraic operations on numerical ones.
     */
    // void operate(const String &expression);

    /** Replace an string in some column(label).
     * The type of the column should be string. This function is a shortcut
     * of the more genereal function operate
     */
    // void replace(const MDLabel label, const String &oldStr, const String &newStr);

    /** Randomize a metadata.
     * MDin is input and the "randomized"
     * result will be in the "calling" Metadata.
    */
    // void randomize(const MetaData &MDin);

    /**Remove duplicate entries for attribute in label
     */
    // void removeDuplicates(MetaData &MDin, MDLabel label=MDL_UNDEFINED);

    /**Remove rows with MDL_ENABLED = -1 if this label is present
     */
    // void removeDisabled();

    /*
    * Sort a Metadata by a label.
    * Sort the content of MDin comparing
    * the label supplied, the result will
    * be in the "calling" MetaData.
    * Limit fixes the maximum number of returned rows
    * Offset skips the first N rows
    */
    /*void sort(MetaData &MDin,
              const MDLabel sortLabel,
              bool asc=true,
              int limit=-1,
              int offset=0);*/


    /*
    * Sort a Metadata by a label.
    * Sort the content of MDin comparing
    * the label supplied, the result will
    * be in the "calling" MetaData.
    * If the input label is a vector field,
    * you may supply label:col, to sort by that column,
    * e.g., NMADisplacements:0
    * Limit fixes the maximum number of returned rows
    * Offset skips the first N rows
    *
    */
    // void sort(MetaData &MDin, const String &sortLabel, bool asc=true, int limit=-1, int offset=0);

    /** Split Metadata in several Metadatas.
     * The Metadata will be divided in 'n'
     * almost equally parts and the result will
     * be a vector of Metadatas. The "calling"
     * Metadata will no be modified.
     * @code
     *   // Divide the images metadata in 10 metadatas.
     *   std::vector<MetaData> imagesGroups;
     *   imageMD.split(10, imagesGroups);
     * @endcode
     */
    //void split(size_t n, std::vector<MetaData> &results,
    //           const MDLabel sortLabel=MDL_OBJID);

    /** Take a part from MetaData.
     * This function is equivallent to divide
     * the input MetaData in n parts and take one.
     * The result will be in "calling" MetaData.
     */
    //void selectSplitPart(const MetaData &mdIn,
    //                     size_t n, size_t part,
    //                     const MDLabel sortLabel=MDL_OBJID);

    /** Select random subset */
    // void selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel=MDL_OBJID);

    /** Select some part from Metadata.
     * Select elements from input Metadata
     * at some starting position
     * if the numberOfObjects is -1, all objects
     * will be returned from startPosition to the end.
    */
    //void selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
    //                const MDLabel sortLabel=MDL_OBJID);

    /** Makes filenames with absolute paths
    *
    */
    // void makeAbsPath(const MDLabel label=MDL_IMAGE);

    /** @} */
    friend struct MDBaseRowIterator<false>;
    friend struct MDBaseRowIterator<true>;

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
        bool operator==(const rowIterator& other) { return other.impl == this->impl; }
        bool operator!=(const rowIterator& other) { return !(*this == other); }
        typename choose<IsConst, const MDRow&, MDRow&>::type operator*() const { return **impl; }
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
        bool operator==(const idIterator& other) { return other.impl == this->impl; }
        bool operator!=(const idIterator& other) { return !(*this == other); }
        size_t operator*() const { return **impl; }
    };

    using id_iterator = idIterator<false>;
    using id_const_iterator = idIterator<true>;

    template <bool IsConst>
    struct IdIteratorProxy {
        typename choose<IsConst, const MetaData&, MetaData&>::type _md;

        IdIteratorProxy(typename choose<IsConst, const MetaData&, MetaData&>::type md) : _md(md) { }
        typename choose<IsConst, MetaData::id_const_iterator, MetaData::id_iterator>::type begin() { return _md.id_begin(); };
        typename choose<IsConst, MetaData::id_const_iterator, MetaData::id_iterator>::type end() { return _md.id_end(); };
    };

    virtual id_iterator id_begin() = 0;
    virtual id_iterator id_end() = 0;

    virtual id_const_iterator id_begin() const = 0;
    virtual id_const_iterator id_end() const = 0;

    virtual IdIteratorProxy<false> ids() { return IdIteratorProxy<false>(*this); };
    virtual IdIteratorProxy<true> ids() const { return IdIteratorProxy<true>(*this); };

    /** Expand Metadata with metadata pointed by label
     * Given a metadata md1, with a column containing the name of another column metdata file mdxx
     * add the columns in mdxx to md1
     */
    // void fillExpand(MDLabel label);

    /** Fill column with constant value
     */
    // void fillConstant(MDLabel label, const String &value);

    /** Fill column with random value
     * mode should be: uniform, gaussian or student
     * op1, op2 and op2 are interpreted for each mode:
     * uniform: op1 and op2 are the limits of the interval
     * gaussian: op1 and op2 are mean and std
     * student: same as gaussian and use op3
     */
    // void fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3=0.);

    /** Fill lineal, starting at some value and with some step */
    // void fillLinear(MDLabel label, double initial, double step);

    /** Copy all values from one column to another.
     * Source column should exist
     */
    // void copyColumn(MDLabel labelDest, MDLabel labelSrc);

    /** Same as previous, but copy to another metadata */
    // void copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc);

    /** Rename column.
     *
     */
    // void renameColumn(MDLabel oldLabel, MDLabel newLabel);

    /** Rename several columns. This is an expensive operations so if several
     * columns need to be changed do it using this function instead one by one
     *
     */
    // void renameColumn(const std::vector<MDLabel> &oldLabel,
    //         const std::vector<MDLabel> &newLabel);

    // void metadataToVec(std::vector<MDRow> &vd);

    // void vecToMetadata(const std::vector<MDRow> &rowMetadata);

    /** 'is equal to' (equality).*/
    // bool operator==(const MetaData& op) const;
}
;//class MetaData

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
