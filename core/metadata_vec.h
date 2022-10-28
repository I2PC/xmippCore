/***************************************************************************
 *
 * Authors:     Jan Horacek (xhorace4@fi.muni.cz)
 *
 * Institute of Computer Science MUNI
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

#ifndef VEC_METADATA_H
#define VEC_METADATA_H

#include <memory>
#include <exception>
#include <unordered_map>
#include "memory.h"
#include "metadata_base.h"
#include "metadata_base_it.h"
#include "utils/memory_utils.h"

using MetaDataVecRow = std::vector<MDObject>;


/** Fast in-memory storage of metadata in std::vector.
 * MetaDataVec implements MetaData API (see metadata_base.h).
 *
 * ### Notes
 *  1. It's fast to iterate either over rows or over ids.
 */
class MetaDataVec: public MetaData {
protected:
    std::vector<MetaDataVecRow> _rows;
    std::array<int, MDL_LAST_LABEL> _label_to_col; // -1 = no mapping
    std::vector<MDLabel> _col_to_label;
    size_t _no_columns = 0;
    std::unordered_map<size_t, size_t> _id_to_index;
    size_t _next_id = 1;

    /** Init, do some initializations tasks, used in constructors
     * @ingroup MetaDataConstructors
     */
    void init(const std::vector<MDLabel>& labelsVector);

    /** clear data and table structure */
    void _clear(bool onlyData=false);

    /**
     * Get a vector of (empty) objects for each active label
     */
    std::vector<MDObject> getObjectsForActiveLabels() const;

    int _labelIndex(MDLabel label) const;
    const MDObject& _getObject(size_t i, MDLabel label) const;
    MDObject& _getObject(size_t i, MDLabel label);
    const MDObject& _getObject(const MetaDataVecRow&, MDLabel) const;
    MDObject& _getObject(MetaDataVecRow&, MDLabel) const;
    int _rowIndex(size_t id) const;
    size_t _rowIndexSafe(size_t id) const;

    void _setRow(const MDRow &row, size_t index);
    bool _match(const MetaDataVecRow&, const MDQuery&) const;
    size_t getRowId(const MetaDataVecRow&) const;

    // Expand row to fit label in.
    void _expand(MetaDataVecRow&, const MDLabel);
    void _expand(MetaDataVecRow&, size_t labeli);

    void _parseObjects(std::istream &is, std::vector<MDObject*> & columnValues,
                       const std::vector<MDLabel> *desiredLabels, bool firstTime) override;

    void _recalc_id_to_index();

    bool _contains(const std::vector<MetaDataVecRow>&, const MetaDataVecRow&) const;
    bool _rowsEq(const MetaDataVecRow& a, const MetaDataVecRow& b) const;

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
    MetaDataVec();
    MetaDataVec(const std::vector<MDLabel> &labelsVector);

    /** From File Constructor.
     *
     * The MetaData is created and data is read from provided FileName. Optionally, a vector
     * of labels can be provided to read just those required labels
     */
    MetaDataVec(const FileName &fileName, const std::vector<MDLabel> &desiredLabels);
    MetaDataVec(const FileName &fileName);

    /** Copy constructor
     *
     * Created a new metadata by copying all data from an existing MetaData object.
     */
    MetaDataVec(const MetaData &md);
    MetaDataVec(const MetaDataVec &md) = default;

    virtual ~MetaDataVec() {}

    /**Clear all data
     */
    void clear() override;
    /** @} */

    /** @name Getters and setters
     * @{
     */
    std::vector<MDLabel> getActiveLabels() const override;

    /** Export medatada to xml file.
     *
     */
    void writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const override;

    /** Write metadata in text file as plain data without header.
     *
     */
    void writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const override;

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */

    size_t addRow(const MDRow &row) override;
    void addRows(const std::vector<MDRowVec> &rows);


    int getMaxStringLength(const MDLabel thisLabel) const override;

    /** Set the value of all objects in an specified column (both value and column are specified in mdValueIn)
    */
    bool setValueCol(const MDObject &mdValueIn) override;

    template<class T>
    bool setValueCol(const MDLabel label, const T &valueIn) {
        return MetaData::setValueCol(label, valueIn);
    }

    /** This functions are using MDObject for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */
    bool setValue(const MDObject &mdValueIn, size_t id);

    template<class T>
    bool setValue(const MDLabel label, const T &valueIn, size_t id) {
        return MetaData::setValue(label, valueIn, id);
    }

    bool getValue(MDObject &mdValueOut, size_t id) const override;
    MDObject &getValue(MDLabel, size_t id);
    const MDObject &getValue(MDLabel, size_t id) const;

    template<class T>
    bool getValue(const MDLabel label, T &valueOut, size_t id) const {
        return MetaData::getValue(label, valueOut, id);
    }

    template<class T>
    T getValue(const MDLabel label, size_t id) {
        return this->getValue(label, id).getValue2(T());
    }

    template<class T>
    const T getValue(const MDLabel label, size_t id) const {
        return this->getValue(label, id).getValue2(T());
    }

    std::unique_ptr<MDRow> getRow(size_t id) override;
    std::unique_ptr<const MDRow> getRow(size_t id) const override;

    MDRowVec getRowVec(size_t id);
    const MDRowVec getRowVec(size_t id) const;

    void getRow(MDRowVec &row, size_t id); // FIXME: deprecated, use getRow above

    bool getRowValues(size_t id, std::vector<MDObject> &values) const override;
    size_t getRowId(size_t i) const;
    void getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const override;
    void setColumnValues(const std::vector<MDObject> &valuesIn) override;

    template<class T>
    void setColumnValues(const MDLabel label, const std::vector<T> &valuesIn) {
        return MetaData::setColumnValues(label, valuesIn);
    }

    bool setRow(const MDRow &row, size_t id);

    template<class T>
    std::vector<T> getColumnValues(const MDLabel label) const {
        return MetaData::getColumnValues<T>(label);
    }

    template<class T>
    void getColumnValues(const MDLabel label, std::vector<T> &valuesOut) const {
        return MetaData::getColumnValues(label, valuesOut);
    }

    /**Check whether the metadata is empty.
     */
    bool isEmpty() const override;

    /**Number of objects contained in the metadata.
     */
    size_t size() const override;

    /** Check whether a label is contained in metadata.
     */
    bool containsLabel(const MDLabel label) const override;

    size_t firstRowId() const override;
    size_t firstObject(const MDQuery&) const override;
    size_t lastRowId() const override;

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */

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
     * A.importObjects(B);
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
    int removeObjects() override;
    int removeObjects(const MDQuery&) override;

    /** @name Search operations
     * @{
     */

    /** Find all objects that match a query.
     * if called without query, all objects are returned
     * if limit is provided only return a maximun of 'limit'
     */
    void findObjects(std::vector<size_t> &objectsOut, const MDQuery &query) const override;
    void findObjects(std::vector<size_t> &objectsOut, int limit = -1) const override;

    size_t countObjects(const MDQuery&) const override;
    bool containsObject(size_t objectId) const override;
    bool containsObject(const MDQuery&) const override;

    /** Find if the object with this id is present in the metadata
     */
    bool containsObject(size_t objectId);

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
    void write(const FileName &outFile, WriteModeMetaData mode=MD_OVERWRITE) const;


    /** Write metadata to out stream
     */
    void write(std::ostream &os, const String & blockName="",WriteModeMetaData mode=MD_OVERWRITE) const;

    /** Read metadata from xml file
     *
     */
    void readXML(const FileName &inFile,
                 const std::vector<MDLabel> *desiredLabels = nullptr,
                 const String & blockRegExp=DEFAULT_BLOCK_NAME,
                 bool decomposeStack=true);

    /** Read data from file. Guess the blockname from the filename
     * @code
     * inFilename="first@md1.doc" -> filename = md1.doc, blockname = first
     * @endcode
     */
    void read(const FileName &inFile, const std::vector<MDLabel> *desiredLabels = nullptr, bool decomposeStack = true) override;
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


    /** Returns Max and Min values from a column in metadata
     * These functions can only be used for labels of type double
     */
    double getColumnMax(MDLabel column);

    double getColumnMin(MDLabel column);


    /** Replace an string in some column(label).
     * The type of the column should be string. This function is a shortcut
     * of the more genereal function operate
     */
    void replace(const MDLabel label, const String &oldStr, const String &newStr);

    /** Randomize a metadata.
     * MDin is input and the "randomized"
     * result will be in the "calling" Metadata.
    */
    void randomize(const MetaData &MDin);

    /**Remove duplicate entries for attribute in label
     */
    void removeDuplicates(MetaData &MDin, MDLabel label=MDL_UNDEFINED);

    /*
    * Sort a Metadata by a label.
    * Sort the content of MDin comparing
    * the label supplied, the result will
    * be in the "calling" MetaData.
    * Limit fixes the maximum number of returned rows
    * Offset skips the first N rows
    */
    void sort(const MetaDataVec &MDin,
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
    void sort(MetaDataVec &MDin, const String &sortLabel, bool asc=true, int limit=-1, int offset=0);

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
    void split(size_t n, std::vector<MetaDataVec> &results,
               const MDLabel sortLabel=MDL_OBJID) const;

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

    /** @name Iterators
     * @{
     */

    template <bool IsConst>
    struct MDVecRowIterator : public MDBaseRowIterator<IsConst> {
    private:
        typename TypeHelpers::choose<IsConst, const MetaDataVec&, MetaDataVec&>::type _mdv;
        size_t _i;
        using RowType = typename TypeHelpers::choose<IsConst, const MDRowVec, MDRowVec>::type;
        std::unique_ptr<RowType> _row;

    public:
        MDVecRowIterator(typename TypeHelpers::choose<IsConst, const MetaDataVec&, MetaDataVec&>::type &mdv, size_t i)
            : _mdv(mdv), _i(i) {
                if (_i < _mdv.size())
                    _row.reset(new RowType(mdv._rows.at(i), i, mdv._label_to_col, mdv._col_to_label, mdv._no_columns));
                else
                    _row = nullptr;
            }

        std::unique_ptr<MDBaseRowIterator<IsConst>> clone() override {
            return memoryUtils::make_unique<MDVecRowIterator<IsConst>>(_mdv, _i);
        }

        void increment() override {
            _i++;
            if (_i < _mdv.size())
                _row.reset(new RowType(_mdv._rows.at(_i), _i, _mdv._label_to_col, _mdv._col_to_label, _mdv._no_columns));
            else
                _row = nullptr;
        }

        bool operator==(const MDBaseRowIterator<IsConst>& other) const override {
            const MDVecRowIterator* vri = dynamic_cast<const MDVecRowIterator<IsConst>*>(&other);
            if (vri != nullptr)
                return _i == vri->_i;
            return false;
        }

        typename TypeHelpers::choose<IsConst, const MDRowVec&, MDRowVec&>::type operator*() override { return *_row; }
    };

    iterator begin() override {
        return {memoryUtils::make_unique<MDVecRowIterator<false>>(*this, 0)};
    }
    iterator end() override {
        return {memoryUtils::make_unique<MDVecRowIterator<false>>(*this, this->size())};
    }

    const_iterator begin() const override {
        return {memoryUtils::make_unique<MDVecRowIterator<true>>(*this, 0)};
    }
    const_iterator end() const override {
        return {memoryUtils::make_unique<MDVecRowIterator<true>>(*this, this->size())};
    }


    template <bool IsConst>
    struct MDVecIdIterator : public MDBaseIdIterator<IsConst> {
    private:
        const MetaDataVec& _mdv;
        size_t _i;

    public:
        MDVecIdIterator(const MetaDataVec& mdv, size_t i)
            : _mdv(mdv), _i(i) {}

        bool operator==(const MDBaseIdIterator<IsConst>& other) const override {
            const MDVecIdIterator<IsConst>* dri = dynamic_cast<const MDVecIdIterator<IsConst>*>(&other);
            if (dri != nullptr)
                return this->_i == dri->_i;
            return false;
        }

        size_t operator*() override { return _mdv.getRowId(_i); }

        void increment() override { this->_i++; }

        std::unique_ptr<MDBaseIdIterator<IsConst>> clone() override {
            return memoryUtils::make_unique<MDVecIdIterator<IsConst>>(_mdv, _i);
        }
    };

    id_iterator id_begin() override {
        return {memoryUtils::make_unique<MDVecIdIterator<false>>(*this, 0)};
    }

    id_iterator id_end() override {
        return {memoryUtils::make_unique<MDVecIdIterator<false>>(*this, this->size())};
    }

    id_const_iterator id_begin() const override {
        return {memoryUtils::make_unique<MDVecIdIterator<true>>(*this, 0)};
    }

    id_const_iterator id_end() const override {
        return {memoryUtils::make_unique<MDVecIdIterator<true>>(*this, this->size())};
    }

    /** @} */

    void fillConstant(MDLabel label, const String &value) override;
    void fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3=0.) override;
    void fillLinear(MDLabel label, double initial, double step) override;

    void copyColumn(MDLabel labelDest, MDLabel labelSrc) override;
    void copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc) override;

    void renameColumn(MDLabel oldLabel, MDLabel newLabel) override;
    void renameColumn(const std::vector<MDLabel> &oldLabel,
            const std::vector<MDLabel> &newLabel) override;

    bool operator==(const MetaDataVec& op) const;
};//class MetaDataVec

std::ostream& operator<<(std::ostream& o, const MetaDataVec& md);

#endif
