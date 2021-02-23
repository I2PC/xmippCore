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
#include "metadata_base.h"
#include "metadata_base_it.h"
#include "metadata_row_vec.h"

class NotImplemented : public std::logic_error
{
public:
    NotImplemented() : std::logic_error("Function not yet implemented") { };
};


using MetaDataVecRow = std::vector<MDObject>;

/** Class to manage data files.
 *
 * The MetaData class manages all procedures related to
 * metadata. MetaData is intended to group together old
 * Xmipp specific files like Docfiles, Selfiles, etc..
 *
 */
class MetaDataVec: public MetaData {
protected:
    std::vector<MetaDataVecRow> _rows;
    std::array<int, MDL_LAST_LABEL> _label_to_col;
    std::vector<MDLabel> _col_to_label;
    size_t _no_columns = 0;

    /** Init, do some initializations tasks, used in constructors
     * @ingroup MetaDataConstructors
     */
    void init(const std::vector<MDLabel>& labelsVector);

    /** Copy info variables from another metadata
     * @ingroup MetaDataConstructors
     */
    void copyInfo(const MetaData &md);

    /** Copy all data from another metadata
     * @ingroup MetaDataConstructors
     */
    void copyMetadata(const MetaData &md, bool copyObjects = true);

    /** This have the same logic of the public one,
     * but doesn't perform any range(which implies do a size()) checks.
     */
    void _selectSplitPart(const MetaData &mdIn,
                          int n, int part, size_t mdSize,
                          const MDLabel sortLabel);

    /** clear data and table structure */
    void _clear(bool onlyData=false);

    /** Some private reading functions */
    void _readColumns(std::istream& is, std::vector<MDObject*> & columnValues,
                      const std::vector<MDLabel>* desiredLabels = nullptr);
    void _readRows(std::istream& is, std::vector<MDObject*> & columnValues, bool useCommentAsImage);
    /** This function will be used to parse the rows data in START format
     * @param[out] columnValues MDRow with values to fill in
     * @param pchStart pointer to the position of '_loop' in memory
     * @param pEnd  pointer to the position of the next '_data' in memory
     * @param maxRows if this number if greater than 0, only this number of rows will be parsed.
     */
    void _readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues, const std::vector<MDLabel> *desiredLabels);
    void _readRowFormat(std::istream& is);

    /**
     * Get a vector of (empty) objects for each active label
     */
    std::vector<MDObject> getObjectsForActiveLabels() const;


    /** This two variables will be used to read the metadata information (labels and size)
     * or maybe a few rows only
     */
    size_t _maxRows, _parsedLines;

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
    MetaDataVec(const MetaDataVec &md);

    /** Assignment operator
     *
     * Copies MetaData from an existing MetaData object.
     */
    MetaDataVec& operator=(const MetaData &md);
    MetaDataVec& operator=(const MetaDataVec &md);

    /** Destructor
     *
     * Frees all used memory and destroys object.
     */
    ~MetaDataVec();

    /**Clear all data
     */
    void clear() override;
    /** @} */

    /** @name Getters and setters
     * @{
     */

    bool nextBlock(mdBuffer &buffer, mdBlock &block);

    /** Export medatada to xml file.
     *
     */
    void writeXML(const FileName fn, const FileName blockname, WriteModeMetaData mode) const override;

    /** Write metadata in text file as plain data without header.
     *
     */
    void writeText(const FileName fn,  const std::vector<MDLabel>* desiredLabels) const override;

    void _parseObjects(std::istream &is, std::vector<MDObject*> & columnValues, const std::vector<MDLabel> *desiredLabels, bool firstTime);

    /* Helper function to parse an MDObject and set its value.
     * The parsing will be from an input stream(istream)
     * and if parsing fails, an error will be raised
     */
    void _parseObject(std::istream &is, MDObject &object, size_t id = BAD_OBJID);

    /** Get Metadata labels for the block defined by start
     * and end loop pointers. Return pointer to newline after last label
     */
    void _readColumnsStar(mdBlock &block,
                          std::vector<MDObject*> & columnValues,
                          const std::vector<MDLabel>* desiredLabels,
                          bool addColumns = true,
                          size_t id = BAD_OBJID);

    /** @} */

    /** @name MetaData Manipulation
     * @{
     */

    size_t addRow(const MDRowVec &row);

    int getMaxStringLength(const MDLabel thisLabel) const override;

    /** Set the value of all objects in an specified column (both value and column are specified in mdValueIn)
    */
    bool setValueCol(const MDObject &mdValueIn) override;

    //private:
    /** This functions are using MDObject for set real values
     * there is an explicit function signature
     * foreach type supported in Metadata.
     * This is done for some type checking of Metadata labels
     * and values
     */
    bool setValue(const MDObject &mdValueIn, size_t id);
    bool getValue(MDObject &mdValueOut, size_t id) const;

    bool getRow(MDRowVec &row, size_t id) const;
    bool getRowValues(size_t id, std::vector<MDObject> &values) const override;
    void getColumnValues(const MDLabel label, std::vector<MDObject> &valuesOut) const override;
    void setColumnValues(const std::vector<MDObject> &valuesIn) override;

    /** Set label values from string representation.
     */
    bool setValueFromStr(const MDLabel label, const String &value, size_t id);

    /** Get string representation from label value.
     */
    bool getStrFromValue(const MDLabel label, String &strOut, size_t id) const;

    /**Check whether the metadata is empty.
     */
    bool isEmpty() const override;

    /**Number of objects contained in the metadata.
     */
    size_t size() const override;

    /** Check whether a label is contained in metadata.
     */
    bool containsLabel(const MDLabel label) const override;

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

    /** Remove all the labels from the metadata but the
     * ones given in labels vector.
     */
    bool keepLabels(const std::vector<MDLabel> &labels) override;

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

    /** Add item id.
     * From 1 to last.
     */
    void addItemId();

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
    void _writeRows(std::ostream &os) const;

    /** Write metadata to disk.
     * This will write the metadata content to disk.
     */
    void writeStar(const FileName &outFile,const String & blockName="", WriteModeMetaData mode=MD_OVERWRITE) const;
    /** Write metadata to disk. Guess blockname from filename
     * @code
     * outFilename="first@md1.doc" -> filename = md1.doc, blockname = first
     * @endcode
     */
    void write(const FileName &outFile, WriteModeMetaData mode=MD_OVERWRITE) const;


    /** Write metadata to out stream
     */
    void write(std::ostream &os, const String & blockName="",WriteModeMetaData mode=MD_OVERWRITE) const;
    void print() const;

    /** Append data lines to file.
     * This function can be used to add new data to
     * an existing metadata. Now should be used with
     * files with only one metadata, maybe can be extended later.
     * For now it will not check any compatibility beetween the
     * existent metadata and the new data to append.
     */
    void append(const FileName &outFile) const;

    /** Check if block exists in metadata file
     * input full parh block@filename
     * return false if metadata block does not exits
     */
    bool existsBlock(const FileName &_inFile);

    /** Read data from file.
     */
    void readStar(const FileName &inFile,
                  const std::vector<MDLabel> *desiredLabels = nullptr,
                  const String & blockName=DEFAULT_BLOCK_NAME,
                  bool decomposeStack=true);
    /** Read metadata from xml file
     *
     */
    void readXML(const FileName &inFile,
                 const std::vector<MDLabel> *desiredLabels= nullptr,
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
    void randomize(const MetaData &MDin);

    /**Remove duplicate entries for attribute in label
     */
    void removeDuplicates(MetaData &MDin, MDLabel label=MDL_UNDEFINED);

    /**Remove rows with MDL_ENABLED = -1 if this label is present
     */
    void removeDisabled();

    /*
    * Sort a Metadata by a label.
    * Sort the content of MDin comparing
    * the label supplied, the result will
    * be in the "calling" MetaData.
    * Limit fixes the maximum number of returned rows
    * Offset skips the first N rows
    */
    void sort(MetaData &MDin,
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
    void sort(MetaData &MDin, const String &sortLabel, bool asc=true, int limit=-1, int offset=0);

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
    void split(size_t n, std::vector<MetaData> &results,
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
    void selectRandomSubset(const MetaData &mdIn, size_t numberOfObjects, const MDLabel sortLabel=MDL_OBJID);

    /** Select some part from Metadata.
     * Select elements from input Metadata
     * at some starting position
     * if the numberOfObjects is -1, all objects
     * will be returned from startPosition to the end.
    */
    void selectPart(const MetaData &mdIn, size_t startPosition, size_t numberOfObjects,
                    const MDLabel sortLabel=MDL_OBJID);

    /** Makes filenames with absolute paths
    *
    */
    void makeAbsPath(const MDLabel label=MDL_IMAGE);

    /** @} */
    template <bool IsConst>
    struct MDVecRowIterator : public MDBaseRowIterator<IsConst> {
    private:
        typename choose<IsConst, const MetaDataVec&, MetaDataVec&>::type _mdv;
        size_t _i;
        typename choose<IsConst, MDRowVecConst, MDRowVec>::type _row;

    public:
        MDVecRowIterator(typename choose<IsConst, const MetaDataVec&, MetaDataVec&>::type &mdv, size_t i)
            : _mdv(mdv), _i(i), _row(mdv._rows.at(i), i, mdv._label_to_col) {}

        // TODO: use std::make_unique when ported to C++14
        std::unique_ptr<MDBaseRowIterator<IsConst>> clone() override {
            return std::unique_ptr<MDVecRowIterator>(new MDVecRowIterator(_mdv, _i));
        }

        void increment() override {
            _i++;
            _row = typename choose<IsConst, MDRowVecConst, MDRowVec>::type(_mdv._rows.at(_i), _i, _mdv._label_to_col);
        }

        bool operator==(const MDBaseRowIterator<IsConst>& other) const override {
            const MDVecRowIterator* vri = dynamic_cast<const MDVecRowIterator<IsConst>*>(&other);
            if (vri != nullptr)
                return _i == vri->_i;
            return false;
        }

        MDRow& operator*() override { return _row; }
    };

    // TODO: use std::make_unique when ported to C++14
    iterator begin() override {
        return rowIterator<false>(std::unique_ptr<MDVecRowIterator<false>>(new MDVecRowIterator<false>(*this, 0)));
    }
    iterator end() override {
        return rowIterator<false>(std::unique_ptr<MDVecRowIterator<false>>(new MDVecRowIterator<false>(*this, this->size())));
    }

    const_iterator begin() const override {
        return rowIterator<true>(std::unique_ptr<MDVecRowIterator<true>>(new MDVecRowIterator<true>(*this, 0)));
    }
    const_iterator end() const override {
        return rowIterator<true>(std::unique_ptr<MDVecRowIterator<true>>(new MDVecRowIterator<true>(*this, this->size())));
    }


    // TODO
    id_iterator id_begin() override {
    }

    id_iterator id_end() override {
    }

    id_const_iterator id_begin() const override {
    }

    id_const_iterator id_end() const override {
    }

    /** Expand Metadata with metadata pointed by label
     * Given a metadata md1, with a column containing the name of another column metdata file mdxx
     * add the columns in mdxx to md1
     */
    void fillExpand(MDLabel label);

    /** Fill column with constant value
     */
    void fillConstant(MDLabel label, const String &value);

    /** Fill column with random value
     * mode should be: uniform, gaussian or student
     * op1, op2 and op2 are interpreted for each mode:
     * uniform: op1 and op2 are the limits of the interval
     * gaussian: op1 and op2 are mean and std
     * student: same as gaussian and use op3
     */
    void fillRandom(MDLabel label, const String &mode, double op1, double op2, double op3=0.);

    /** Fill lineal, starting at some value and with some step */
    void fillLinear(MDLabel label, double initial, double step);

    /** Copy all values from one column to another.
     * Source column should exist
     */
    void copyColumn(MDLabel labelDest, MDLabel labelSrc);

    /** Same as previous, but copy to another metadata */
    void copyColumnTo(MetaData& md, MDLabel labelDest, MDLabel labelSrc);

    /** Rename column.
     *
     */
    void renameColumn(MDLabel oldLabel, MDLabel newLabel);

    /** Rename several columns. This is an expensive operations so if several
     * columns need to be changed do it using this function instead one by one
     *
     */
    void renameColumn(const std::vector<MDLabel> &oldLabel,
            const std::vector<MDLabel> &newLabel);

    void metadataToVec(std::vector<MDRow> &vd);

    void vecToMetadata(const std::vector<MDRow> &rowMetadata);

    /** 'is equal to' (equality).*/
    bool operator==(const MetaData& op) const;
}
;//class MetaData

/** print metadata
 *
 */
std::ostream& operator<<(std::ostream& o, const MetaData & mD);

#endif
