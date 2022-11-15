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

#ifndef CORE_METADATA_DB_H
#define CORE_METADATA_DB_H

#include <regex.h>
#include <cmath>
#include "utils/memory_utils.h"
#include "metadata_base.h"
#include "metadata_label.h"
#include "metadata_object.h"
#include "metadata_row_base.h"
#include "metadata_static.h"
#include "metadata_sql.h"
#include "metadata_sql_operations.h"
#include "utils/sql_utils.h"
#include "xmipp_error.h"
#include "xmipp_filename.h"
#include "metadata_writemode.h"


/** Original database implementation of MetaData.
 * MetaData are stored in SQL database.
 * Some database-specific commands are available in this implementation only.
 *
 * ### Notes
 *  1. It's quite fast to iterate over ids.
 *  2. It's slow to iterate over rows.
 *  3. Best practices: <https://github.com/I2PC/xmipp-portal/wiki/MetaData---SQL-best-practices>
 */
class MetaDataDb : public MetaData {
protected:
    /** This variables should only be used by MDSql
     * for handling db status of metadata
     */
    /** The table id to do db operations */
    friend class MDSql;
    MDSql * myMDSql;

    /** What labels have been read from a docfile/metadata file
     * and/or will be stored on a new metadata file when "save" is
     * called
     **/
    std::vector<MDLabel> _activeLabels;

    /** Init, do some initializations tasks, used in constructors
     * @ingroup MetaDataConstructors
     */
    void init(const std::vector<MDLabel> &labelsVector);

    /** Copy all data from another metadata
     * @ingroup MetaDataConstructors
     */
    void copyMetadata(const MetaDataDb &md, bool copyObjects = true);

    /** This have the same logic of the public one,
     * but doesn't perform any range(which implies do a size()) checks.
     */
    void _selectSplitPart(const MetaDataDb &mdIn,
                          int n, int part, size_t mdSize,
                          const MDLabel sortLabel);

    void _selectSplitPart(const MetaDataDb &mdIn,
                          size_t n, size_t part,
                          const MDLabel sortLabel=MDL_OBJID);

    void _selectRandomSubset(const MetaDataDb &mdIn, size_t numberOfObjects, const MDLabel sortLabel=MDL_OBJID);

    void _selectPart(const MetaDataDb &mdIn, size_t startPosition, size_t numberOfObjects,
                     const MDLabel sortLabel=MDL_OBJID);

    /** This function is for generalize the sets operations
     * of unionDistinct, intersection, subtraction
     * which can be expressed in terms of
     * ADD, SUBSTRACT of intersection part
     */
    void _setOperates(const MetaDataDb &mdIn, const MDLabel label, SetOperation operation);
    void _setOperates(const MetaDataDb &mdIn, const std::vector<MDLabel> &labels, SetOperation operation);
    void _setOperates(const MetaDataDb &mdInLeft,
                      const MetaDataDb &mdInRight,
                      const std::vector<MDLabel> &labelsLeft,
                      const std::vector<MDLabel> &labelsRight,
                      SetOperation operation);
    /** This function is for generalize the sets operations
     * in which the output has a single label
     * a vector of labels instead of a single label may be implemented in the future
     */
    void _setOperatesLabel(const MetaDataDb &mdIn, const MDLabel label, SetOperation operation);
    /** clear data and table structure */
    void _clear(bool onlyData=false);

    void _readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues,
                       const std::vector<MDLabel> *desiredLabels) override;

    /** Some private reading functions */
    void _readRowFormat(std::istream& is);

    void _parseObjects(std::istream &is, std::vector<MDObject*> & columnValues,
                       const std::vector<MDLabel> *desiredLabels, bool firstTime) override;

    /**
     * Get a vector of (empty) objects for each active label
     */
    std::vector<MDObject> getObjectsForActiveLabels() const;

    void _importObjectsDb(const MetaDataDb &md, const MDQuery &query, bool doClear=true);
    void _importObjectsGeneral(const MetaData &md, const MDQuery &query, bool doClear=true);

public:
    /** @name Constructors
     *  @{
     */

    /** Empty Constructor.
     *
     * The MetaDataDb is created with no data stored on it. You can fill in it programmatically
     * or by a later reading from a MetaDataDb file or old Xmipp formatted type.
     * if labels vectors is passed this labels are created on metadata
     */
    MetaDataDb();
    MetaDataDb(const std::vector<MDLabel> &labelsVector);
    MetaDataDb(const MetaData &md);

    /** From File Constructor.
     *
     * The MetaData is created and data is read from provided FileName. Optionally, a vector
     * of labels can be provided to read just those required labels
     */
    MetaDataDb(const FileName &fileName, const std::vector<MDLabel> &desiredLabels = {});

    /** Copy constructor
     *
     * Created a new metadata by copying all data from an existing MetaData object.
     */
    MetaDataDb(const MetaDataDb &md);

    /** Assignment operator
     *
     * Copies MetaDataDb from an existing MetaData object.
     */
    MetaDataDb& operator=(const MetaDataDb &md);

    /** Destructor
     *
     * Frees all used memory and destroys object.
     */
    virtual ~MetaDataDb();

    /**Clear all data
     */
    void clear() override;
    /** @} */

    /** @name Getters and setters
     * @{
     */

    MDSql * getDatabase() { return myMDSql; }

    /** Export medatada to xml file.
     *
     */
    void writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const override;

    /** Write metadata in sqlite3 file.
     *
     */
    void writeDB(const FileName fn, const FileName blockname, WriteModeMetaData mode) const;

    /** Write metadata in text file as plain data without header.
     *
     */
    void writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const override;

    /**Get maximum string length of column values.
    */
    int getMaxStringLength( const MDLabel thisLabel) const override;

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */


    /** Set the value of all objects in an specified column (both value and column are specified in mdValueIn)
    */
    bool setValueCol(const MDObject &mdValueIn) override;

    template<class T>
    bool setValueCol(const MDLabel label, const T &valueIn) {
        return MetaData::setValueCol(label, valueIn);
    }

    //private:
    /** This functions are using MDObject for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */
    bool setValue(const MDObject &mdValueIn, size_t id) override;


    template<class T>
    bool setValue(const MDLabel label, const T &valueIn, size_t id) {
        return MetaData::setValue(label, valueIn, id);
    }

    bool getValue(MDObject &mdValueOut, size_t id) const override;
    bool getRowValues(size_t id, std::vector<MDObject> &values) const override;

    template<class T>
    bool getValue(const MDLabel label, T &valueOut, size_t id) const {
        return MetaData::getValue(label, valueOut, id);
    }

    /** Get all values of a column as a vector.
     */
    void getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const override;

    template<class T>
    std::vector<T> getColumnValues(const MDLabel label) const {
        return MetaData::getColumnValues<T>(label);
    }

    template<class T>
    void getColumnValues(const MDLabel label, std::vector<T> &valuesOut) const {
        return MetaData::getColumnValues(label, valuesOut);
    }

    /** Get all values of a column as a vector.
     */
    template<typename T>
    bool getColumnValuesOpt(const MDLabel label, std::vector<T> &values) const;

    /** Get all values of a column as a vector.
     */
    void setColumnValues(const std::vector<MDObject> &valuesIn) override;

    template<class T>
    void setColumnValues(const MDLabel label, const std::vector<T> &valuesIn) {
        return MetaData::setColumnValues(label, valuesIn);
    }

    /** Get all values of an MetaData row of an specified objId*/
    bool bindValue(size_t id) const;

    bool initGetRow(bool addWhereClause) const;
    bool execGetRow(MDRow &row) const;
    void finalizeGetRow(void) const;

    std::unique_ptr<MDRow> getRow(size_t id) override;
    std::unique_ptr<const MDRow> getRow(size_t id) const override;

    MDRowSql getRowSql(size_t id);
    const MDRowSql getRowSql(size_t id) const;

    bool getRow(MDRowSql &row, size_t id) const; // FIXME: deprecated, use getRow above

    bool getAllRows(std::vector<MDRowSql> &rows) const;
    bool getRow2(MDRow &row, size_t id) const;

    /** Copy all the values in the input row in the current metadata*/
    bool setRow(const MDRow &row, size_t id);

    /** Add a new Row and set values, return the objId of newly added object */
    bool initAddRow(const MDRow &row);
    bool execAddRow(const MDRow &row);
    void finalizeAddRow(void);
    size_t addRow(const MDRow &row) override;
    void addRowOpt(const MDRowSql &row);
    void addRows(const std::vector<MDRowSql> &rows);
    void addMissingLabels(const MDRow &row);
    size_t addRow2(const MDRow &row);

    /**Number of objects contained in the metadata.
     */
    size_t size() const override;

    bool containsLabel(const MDLabel label) const override {
        return vectorContainsLabel(this->_activeLabels, label);
    }

    std::vector<MDLabel> getActiveLabels() const override {
        return this->_activeLabels;
    }

    /** Add a new label to the metadata.
     * By default the label is added at the end,
     * if the position is specified and is between 0 and n-1
     * the new label is inserted at that position.
     */
    bool addLabel(const MDLabel label, int pos = -1) override;

    /** Remove a label from the metadata.
     * The data is still in the table. If you want to remove the data,
     * make a copy of the MetaData.
     */
    bool removeLabel(const MDLabel label) override;

    /** Adds a new, empty object to the objects map. If objectId == -1
     * the new ID will be that for the last object inserted + 1, else
     * the given objectId is used. If there is already an object whose
     * objectId == input objectId, just removes it and creates an empty
     * one
     */
    size_t addObject() override;

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
    void importObject(const MetaData &md, const size_t id, bool doClear=true) override;
    void importObjects(const MetaData &md, const std::vector<size_t> &objectsToAdd, bool doClear=true) override;
    void importObjects(const MetaData &md, const MDQuery &query, bool doClear=true) override;

    /** Remove the object with this id.
    * Returns true if the object was removed or false if
    * the object did not exist
    */
    bool removeObject(size_t id) override;

    /** Removes the collection of objects of given vector id's
     * NOTE: The iterator will point to the first object after any of these
     * operations
     */
    void removeObjects(const std::vector<size_t> &toRemove) override;

    /** Removes objects from metadata.
     * return the number of deleted rows
     * if not query, all objectes are removed
     * Queries can be used in the same way
     * as in the importObjects function
     */
    int removeObjects(const MDQuery &query) override;
    int removeObjects() override;

    /** Add and remove indexes for fast search
     * in other labels, but insert are more expensive
     */
    void addIndex(MDLabel label) const;
    void addIndex(const std::vector<MDLabel> &desiredLabels) const;
    void removeIndex(MDLabel label);
    void removeIndex(const std::vector<MDLabel> &desiredLabels);

    /** Add item id.
     * From 1 to last.
     */
    void addItemId();

    /** Remove item id.*/
    void removeItemId();

    /** @} */

    /** @name Iteration functions
     * @{
     */

    size_t firstRowId() const override;
    size_t firstObject(const MDQuery&) const override;
    size_t lastRowId() const override;

    /** @name Search operations
     * @{
     */

    /** Find all objects that match a query.
     * if called without query, all objects are returned
     * if limit is provided only return a maximun of 'limit'
     */
    void findObjects(std::vector<size_t> &objectsOut, const MDQuery &query) const override;
    void findObjects(std::vector<size_t> &objectsOut, int limit = -1) const override;

    size_t countObjects(const MDQuery &query) const override;
    bool containsObject(size_t objectId) const override;
    bool containsObject(const MDQuery &query) const override;

    /** @} */

    /** @name I/O functions
     * @{
     */

    /** Write rows data to disk. */
    void _writeRows(std::ostream &os) const override;

    /** Write metadata to disk. Guess blockname from filename
     * @code
     * outFilename="first@md1.doc" -> filename = md1.doc, blockname = first
     * @endcode
     */
    void write(const FileName &outFile, WriteModeMetaData mode=MD_OVERWRITE) const override;
    void write(std::ostream &os, const String & blockName="",WriteModeMetaData mode=MD_OVERWRITE) const override;

    /** Check if block exists in metadata file
     * input full parh block@filename
     * return false if metadata block does not exits
     */
    bool existsBlock(const FileName &_inFile);

    /** Read metadata from xml file
     *
     */
    void readXML(const FileName &inFile,
                 const std::vector<MDLabel> *desiredLabels= NULL,
                 const String & blockRegExp=DEFAULT_BLOCK_NAME,
                 bool decomposeStack=true);

    /** Read metadata from sqlite file
     *
     */
    void readDB(const FileName &inFile,
                const std::vector<MDLabel> *desiredLabels= NULL,
                const String & blockRegExp=DEFAULT_BLOCK_NAME,
                bool decomposeStack=true);

    /** Read data from file. Guess the blockname from the filename
     * @code
     * inFilename="first@md1.doc" -> filename = md1.doc, blockname = first
     * @endcode
     */
    void read(const FileName &inFile, const std::vector<MDLabel> *desiredLabels = NULL, bool decomposeStack=true) override;
    /** @} */

    /** Try to read a metadata from plain text with some columns.
     * Labels for each columns should be provided in an string separated by spaces.
     *  Return false if couldn't read
     */
    void readPlain(const FileName &inFile, const String &labelsString, const String &separator = " ");
    /** Same as readPlain, but instead of cleanning data, the
     * readed values will be added. If there are common columns in metadata
     * and the plain text, the lattest will be setted
     */
    void addPlain(const FileName &inFile, const String &labelsString, const String &separator=" ");

    /** @name Set Operations
     * @{
     */

    /** Aggregate metadata objects,
     * result in calling metadata object (except for aggregateSingle)
     * thisLabel label is used for aggregation, second. Valid operations are:
     *
     * MDL_AVG:  The avg function returns the average value of all  operationLabel within a group.
      The result of avg is always a floating point value as long as at there
      is at least one non-NULL input even if all inputs are integers.
       The result of avg is NULL if and only if there are no non-NULL inputs.

      AGGR_COUNT: The count function returns a count of the number of times that operationLabel is in a group.

      AGGR_MAX       The max aggregate function returns the maximum value of all values in the group.

      AGGR_MIN       The min aggregate function returns the minimum  value of all values in the group.

     AGGRL_SUM The total aggregate functions return sum of all values in the group.
     If there are no non-NULL input rows then returns 0.0.


     */

    void aggregate(const MetaDataDb &mdIn, AggregateOperation op,
                   MDLabel aggregateLabel, MDLabel operateLabel, MDLabel resultLabel);
    void aggregate(const MetaDataDb &mdIn, const std::vector<AggregateOperation> &ops,
                   const std::vector<MDLabel> &operateLabels, const std::vector<MDLabel> &resultLabels);
    void aggregateGroupBy(const MetaDataDb &mdIn,
                          AggregateOperation op,
                          const std::vector<MDLabel> &groupByLabels,
                          MDLabel operateLabel,
                          MDLabel resultLabel);
    /** This function performs aggregation operations.
        without grouping. (i.e. absolute maximum of a metadata column)
        for double
     */
    void aggregateSingle(MDObject &mdValueOut, AggregateOperation op,
                         MDLabel aggregateLabel);
    /** This function performs aggregation operations.
        without grouping. (i.e. absolute maximum of a metadata column)
        for int
     */
    void aggregateSingleInt(MDObject &mdValueOut, AggregateOperation op,
                            MDLabel aggregateLabel);
    /** This function performs aggregation operations.
        without grouping. (i.e. absolute maximum of a metadata column)
        for size_t
     */
    void aggregateSingleSizeT(MDObject &mdValueOut, AggregateOperation op,
                              MDLabel aggregateLabel);



    /** Returns Max and Min values from a column in metadata
     * These functions can only be used for labels of type double
     */
    double getColumnMax(MDLabel column);

    double getColumnMin(MDLabel column);


    /** Union of elements in two Metadatas, without duplicating.
     * Result in calling metadata object
     * union is a reserved word so I called this method unionDistinct
     */
    void unionDistinct(const MetaDataDb &mdIn, const MDLabel label=MDL_OBJID);

    /** Union of all elements in two Metadata, duplicating common elements.
     * Result in calling metadata object
     * Repetition are allowed
     */
    void unionAll(const MetaDataDb &mdIn);

    /** Merge of two Metadata.
     * This function reads another metadata and add all columns values.
     * The size of the two Metadatas should be the same. If there are
     * common columns, the values in md2 will be setted.
     */
    void merge(const MetaData &md2);

    /** Intersects two Metadatas.
     * Result in "calling" Metadata
     */
    void intersection(const MetaDataDb &mdIn, const MDLabel label);

    /** Subtract two Metadatas.
     * Result in "calling" metadata
     */
    void subtraction(const MetaDataDb &mdIn, const MDLabel label);

    /** Return only distinct (different) values of column label.
     * Result in "calling" metadata with a single column
     */
    void distinct(MetaDataDb &MDin, MDLabel label);

    /** Join two Metadatas
     * Result in "calling" metadata
     */
    void join1(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const MDLabel label, JoinType type=LEFT);

    /** Join two Metadatas
     * Result in "calling" metadata. join may be done using different labels in each metadata
     */
    void join2(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const MDLabel labelLeft, const MDLabel labelRight , JoinType type=LEFT);

    /** Join two Metadatas
     * Result in "calling" metadata
     */
    void join1(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const std::vector<MDLabel> &labels, JoinType type=LEFT);

    /** Join two Metadatas
     * Result in "calling" metadata. join may be done using different labels in each metadata
     */
    void join2(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight, const std::vector<MDLabel> &labelsLeft, const std::vector<MDLabel> &labelsRight,
    		JoinType type=LEFT);

    /** Join two Metadatas using all common labels (NATURAL_JOIN)
     */
    void joinNatural(const MetaDataDb &mdInLeft, const MetaDataDb &mdInRight);

    /** Basic operations on columns data.
     * Mainly perform replacements on string values and
     * basic algebraic operations on numerical ones.
     */
    void operate(const String &expression);

    /** Replace an string in some column(label).
     * The type of the column should be string. This function is a shortcut
     * of the more genereal function operate
     */
    void replace(const MDLabel label, const String &oldStr, const String &newStr);

    /** Randomize a metadata.
     * MDin is input and the "randomized"
     * result will be in the "calling" Metadata.
    */
    void randomize(const MetaDataDb &MDin);

    /**Remove duplicate entries for attribute in label
     */
    void removeDuplicates(MetaDataDb &MDin, MDLabel label=MDL_UNDEFINED);

    /*
    * Sort a Metadata by a label.
    * Sort the content of MDin comparing
    * the label supplied, the result will
    * be in the "calling" MetaData.
    * Limit fixes the maximum number of returned rows
    * Offset skips the first N rows
    */
    void sort(MetaDataDb &MDin,
              const MDLabel sortLabel,
              bool asc=true,
              int limit=-1,
              int offset=0);


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
    void sort(MetaDataDb &MDin, const String &sortLabel, bool asc=true, int limit=-1, int offset=0);

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
    void split(size_t n, std::vector<MetaDataDb> &results,
               const MDLabel sortLabel=MDL_OBJID);

    /** Take a part from MetaData.
     * This function is equivallent to divide
     * the input MetaData in n parts and take one.
     * The result will be in "calling" MetaData.
     */
    void selectSplitPart(const MetaData &mdIn,
                         size_t n, size_t part,
                         const MDLabel sortLabel=MDL_OBJID);

    /** Select random subset */
    void selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel=MDL_OBJID) override;

    /** Select some part from Metadata.
     * Select elements from input Metadata
     * at some starting position
     * if the numberOfObjects is -1, all objects
     * will be returned from startPosition to the end.
    */
    void selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                    const MDLabel sortLabel=MDL_OBJID) override;

    /** Makes filenames with absolute paths
    *
    */
    void makeAbsPath(const MDLabel label=MDL_IMAGE);

    /** @} */

    template <bool IsConst>
    struct MDDbRowIterator : public MDBaseRowIterator<IsConst> {
    private:
        typename TypeHelpers::choose<IsConst, const MetaDataDb&, MetaDataDb&>::type _mdd;
        typename TypeHelpers::choose<IsConst, MDRowSql, MDRowSql>::type _row;
        std::vector<size_t> _ids;
        size_t _i;
        bool _finalized = false;

    public:
        MDDbRowIterator(typename TypeHelpers::choose<IsConst, const MetaDataDb&, MetaDataDb&>::type &mdd, size_t _i)
            : _mdd(mdd), _i(_i) {
            mdd.myMDSql->selectObjects(_ids);
            if (this->_i >= _ids.size())
                return;

            _mdd.initGetRow(false);

            _mdd.execGetRow(this->_row);
            this->_row.set_id(this->_ids[this->_i]);
            if (this->_i+1 == _ids.size()) {
                _mdd.finalizeGetRow();
                _finalized = true;
            }
        }

        virtual ~MDDbRowIterator() {
            if (!_finalized)
                _mdd.finalizeGetRow();
        }

        std::unique_ptr<MDBaseRowIterator<IsConst>> clone() override {
            return memoryUtils::make_unique<MDDbRowIterator<IsConst>>(_mdd, _i);
        }

        void increment() override {
            if (this->_i >= _ids.size())
                return;

            this->_i++;

            if (this->_i >= _ids.size())
                return;

            this->_mdd.execGetRow(this->_row);
            this->_row.set_id(this->_ids[this->_i]);

            if (this->_i == _ids.size()) {
                _mdd.finalizeGetRow();
                _finalized = true;
            }
        }

        bool operator==(const MDBaseRowIterator<IsConst>& other) const override {
            const MDDbRowIterator<IsConst>* dri = dynamic_cast<const MDDbRowIterator<IsConst>*>(&other);
            if (dri != nullptr)
                return this->_i == dri->_i;
            return false;
        }

        typename TypeHelpers::choose<IsConst, const MDRowSql&, MDRowSql&>::type operator*() override { return _row; }
    };

    iterator begin() override {
        return {memoryUtils::make_unique<MDDbRowIterator<false>>(*this, 0)};
    }
    iterator end() override {
        return {memoryUtils::make_unique<MDDbRowIterator<false>>(*this, this->size())};
    }

    const_iterator begin() const override {
        return {memoryUtils::make_unique<MDDbRowIterator<true>>(*this, 0)};
    }
    const_iterator end() const override {
        return {memoryUtils::make_unique<MDDbRowIterator<true>>(*this, this->size())};
    }


    template <bool IsConst>
    struct MDDbIdIterator : public MDBaseIdIterator<IsConst> {
    private:
        const MetaDataDb& _mdd;
        bool _last;
        const MDQuery* _pQuery;
        std::vector<size_t> _ids;
        size_t _i;

    public:
        MDDbIdIterator(const MetaDataDb& mdd, bool last = false, const MDQuery* pQuery = nullptr)
            : _mdd(mdd), _last(last), _pQuery(pQuery) {
            mdd.myMDSql->selectObjects(_ids, pQuery);
            _i = last ? this->_ids.size() : 0;
        }

        bool operator==(const MDBaseIdIterator<IsConst>& other) const override {
            const MDDbIdIterator<IsConst>* dri = dynamic_cast<const MDDbIdIterator<IsConst>*>(&other);
            if (dri != nullptr)
                return this->_i == dri->_i;
            return false;
        }

        size_t operator*() override { return _ids[_i]; }

        void increment() override { this->_i++; }

        std::unique_ptr<MDBaseIdIterator<IsConst>> clone() override {
            return memoryUtils::make_unique<MDDbIdIterator<IsConst>>(_mdd, _last, _pQuery);
        }
    };

    id_iterator id_begin() override {
        return {memoryUtils::make_unique<MDDbIdIterator<false>>(*this)};
    }

    id_iterator id_end() override {
        return {memoryUtils::make_unique<MDDbIdIterator<false>>(*this, true)};
    }

    id_const_iterator id_begin() const override {
        return {memoryUtils::make_unique<MDDbIdIterator<true>>(*this)};
    }

    id_const_iterator id_end() const override {
        return {memoryUtils::make_unique<MDDbIdIterator<true>>(*this, true)};
    }

    /** Expand Metadata with metadata pointed by label
     * Given a metadata md1, with a column containing the name of another column metdata file mdxx
     * add the columns in mdxx to md1
     */
    void fillExpand(MDLabel label);

    void fillConstant(MDLabel label, const String &value) override;
    void fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3=0.) override;
    void fillLinear(MDLabel label, double initial, double step) override;

    void copyColumn(MDLabel labelDest, MDLabel labelSrc) override;
    void copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc) override;

    void renameColumn(MDLabel oldLabel, MDLabel newLabel) override;
    void renameColumn(const std::vector<MDLabel> &oldLabel,
            const std::vector<MDLabel> &newLabel) override;

    void metadataToVec(std::vector<MDRowSql> &vd);
    void vecToMetadata(const std::vector<MDRow> &rowMetadata);

    /** 'is equal to' (equality).*/
    bool operator==(const MetaDataDb& op) const;
}
;//class MetaData

/** print metadata
 *
 */
std::ostream& operator<<(std::ostream& o, const MetaData & mD);

/** Convert string to write mode metadata enum.
 *
 */
WriteModeMetaData metadataModeConvert (String mode);

template<typename T>
bool MetaDataDb::getColumnValuesOpt(const MDLabel label, std::vector<T> &values) const {
    if (!containsLabel(label))
            return false;
    return sqlUtils::select(label,
            myMDSql->db,
            myMDSql->tableName(myMDSql->tableId),
            values);
}

/** @} */

#endif
