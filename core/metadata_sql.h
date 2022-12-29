/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

/* SQL helpers for MetaDataDb */

#ifndef CORE_METADATASQL_H
#define CORE_METADATASQL_H

#include "xmipp_strings.h"
#include "metadata_label.h"
#include "xmipp_error.h"
#include "metadata_sql_operations.h"
#include "metadata_static.h"
#include "metadata_query.h"

class MDSqlStaticInit;
class MDQuery;
class MetaDataDb;
class MDCache;
class FileName;
class sqlite3_stmt;
class sqlite3_context;
class sqlite3_value;
class sqlite3;

/** @addtogroup MetaData
 * @{
 */

//#include "metadata.h"

/* return number of tables from a metadata file saved as sqlite */
int getBlocksInMetaDataFileDB(const FileName &inFile, StringVector& blockList);

/*support for the REGEXP operator in sqlite*/
void sqlite_regexp(sqlite3_context* context, int argc, sqlite3_value** values);

/** This class will manage SQL database interactions.
 * This class is designed to used inside a MetaData.
 */
class MDSql
{
public:
    static void dumpToFile(const FileName &fileName);
    static void sqlTimeOut(int miliSeconds);

    /**This library will provide common mathematical and string functions in
SQL queries using the operating system libraries or provided
definitions.  It includes the following functions:

Math: acos, asin, atan, atn2, atan2, acosh, asinh, atanh, difference,
degrees, radians, cos, sin, tan, cot, cosh, sinh, tanh, coth, exp,
log, log10, power, sign, sqrt, square, ceil, floor, pi.

String: replicate, charindex, leftstr, rightstr, ltrim, rtrim, trim,
replace, reverse, proper, padl, padr, padc, strfilter.

Aggregate: stdev, variance, mode, median, lower_quartile,
upper_quartile.
     *
     *
     */
    static bool activateMathExtensions(void);

    /* activate regular expressions in sql */
    static bool activateRegExtensions(void);

/** activate UNSAFED Multi-thread mode
 * SQLite support three different threading modes:

Single-thread. In this mode, all mutexes are disabled and SQLite is unsafe to use in more than a single thread at once.

Multi-thread. In this mode, SQLite can be safely used by multiple threads provided that no single database connection is used simultaneously in two or more threads.

Serialized. In serialized mode, SQLite can be safely used by multiple threads with no restriction.

 */
    bool  deactivateThreadMuting(void);
    /** activate SAFED serialized mode mode
*/
    bool  activateThreadMuting(void);

private:
    /** write metadata in sqlite table
     *
     */
    void copyTableToFileDB(const FileName blockname, const FileName &fileName);

    /** read metadata from sqlite table
     *
     */
    void copyTableFromFileDB(const FileName blockname,
                             const FileName filename,
                             const std::vector<MDLabel> *desiredLabels,
                             const size_t maxRows=0
                             );
    /** This will create the table to store the metada objects.
     * Will return false if the mdId table is already present.
     */
    bool createMd();

    /** This function will drop the entire table.
     * For use the metada again, a call to createMd() should be done.
     */
    bool clearMd();

    size_t getObjId();

    /**Add a new row and return the objId(rowId).
     */
    size_t addRow();

    /** Add a new column to a metadata.
     */
    bool addColumn(MDLabel column);

   /** Rename Column
    * SQLite itself does not support it. So some hacking is needed here
    */
    bool renameColumn(const std::vector<MDLabel> &oldLabel, const std::vector<MDLabel> &newlabel);

    /** Insert a new register inserting input columns.
     */
    // T is either "const MDObject*" or "MDObject*"
    template <typename T>
    bool setObjectValues(int id, const std::vector<T> &columnValues, const std::vector<MDLabel> *desiredLabels=NULL);

    /**Set the value of an object in an specified column.
     */
    bool setObjectValue(const int objId, const MDObject &value);

    /**Set the value of all objects in an specified column.
     */
    bool setObjectValue(const MDObject &value);

    /** Get the values of several objects.
     */
    bool getObjectsValues(const std::vector<MDLabel> &labels, std::vector<MDObject> &values);

    /** Get the value of an object.
     */
    bool getObjectValue(const int objId, MDObject  &value);

    /** This function will select some elements from table.
     * The 'limit' is the maximum number of object
     * returned, if is -1, all will be returned
     * Also a query could be specified for selecting objects
     * if no query is provided by default all are returned
     */
    void selectObjects(std::vector<size_t> &objectsOut, const MDQuery *queryPtr = NULL);

    /** return metadata size
     *
     */
    size_t size(void);
    /** This function will delete elements that match the query.
     * If not query is provided, all rows are deleted
     */
    size_t deleteObjects(const MDQuery *queryPtr = NULL);

    /** Copy the objects from a metada to other.
     * return the number of objects copied
     * */
    size_t copyObjects(MDSql * sqlOut,
                       const MDQuery *queryPtr = NULL) const;
    size_t copyObjects(MetaDataDb * mdPtrOut,
                       const MDQuery *queryPtr = NULL) const;

    /** This function performs aggregation operations.
     */
    void aggregateMd(MetaDataDb *mdPtrOut,
                     const std::vector<AggregateOperation> &operations,
                     const std::vector<MDLabel> &operateLabel);

    /** This function performs aggregation operations grouped by several labels.
     */
    void aggregateMdGroupBy(MetaDataDb *mdPtrOut,
                            const AggregateOperation operation,
                            const std::vector<MDLabel> &groupByLabels ,
                            const  MDLabel operateLabel,
                            const  MDLabel resultLabel);

    /** This function performs aggregation operations.
        without grouping. (i.e. absolute maximum of a metadata column)
        for double
     */
    double aggregateSingleDouble(const AggregateOperation operation,
                                 MDLabel operateLabel);


    /** This function performs aggregation operations.
        without grouping. (i.e. absolute maximum of a metadata column)
        for size_t
     */
    size_t aggregateSingleSizeT(const AggregateOperation operation,
                                MDLabel operateLabel);

    /** This function will be used to create o delete an index over a column.
     *Those indexes will improve searchs, but inserts will become expensives
     */
    void indexModify(const std::vector<MDLabel> &columns, bool create=true);

    /** Some iteration methods
     */
    size_t firstRow();
    size_t lastRow();
    size_t nextRow(size_t currentRow);
    size_t previousRow(size_t currentRow);

    int columnMaxLength(MDLabel column);

    /**Functions to implement set operations */
    void setOperate(MetaDataDb *mdPtrOut, const std::vector<MDLabel> &columns, SetOperation operation);
    void setOperate(const MetaDataDb *mdInLeft, const MetaDataDb *mdInRight, const std::vector<MDLabel> &columnsLeft,
    		const std::vector<MDLabel> &columnsRight, SetOperation operation);
    /** Function to dump DB to file */
    bool operate(const String &expression);

    /** If true, all queries to database will be synchronized by mutex */
    // FIXME this should not be necessary, as sqlite is build with multi-thread support. However, multiple access
    // to database from multiple threads causes application crash
    bool beThreadSafe;


    /** Constructor of MDSql
     * Now each MD should have an instance
     * of this class to interact with the DB
     */
    MDSql(MetaDataDb *md);
    ~MDSql();

    static int table_counter;
    static sqlite3 *db;

    static MDSqlStaticInit initialization; //Just for initialization

    ///Just call this function once, at static initialization
    static bool sqlBegin();
    static void sqlEnd();
    static bool sqlBeginTrans();
    static bool sqlCommitTrans();
    /** Return an unique id for each metadata
     * this function should be called once for each
     * metada and the id will be used for operations
     */
    int getUniqueId();

    bool dropTable();
    bool createTable(const std::vector<MDLabel> * labelsVector = NULL, bool withObjID=true);
    bool insertValues(double a, double b);
    bool initializeSelect( bool addWhereObjId, const std::vector<MDLabel> &labels);
    bool initializeInsert(const std::vector<MDLabel> *labels, const std::vector<MDObject*> &values);
    void finalizePreparedStmt(void);
    void prepareStmt(const std::stringstream &ss, sqlite3_stmt *stmt);
    bool execSingleStmt(const std::stringstream &ss);
    bool execSingleStmt(sqlite3_stmt *&stmt, const std::stringstream *ss = NULL);
    size_t execSingleIntStmt(const std::stringstream &ss);
    double execSingleDoubleStmt(const std::stringstream &ss);

    String tableName(const int tableId) const;

    bool 	bindStatement( size_t id);
    int 	bindValue(sqlite3_stmt *stmt, const int position, const MDObject &valueIn);
    void 	extractValue(sqlite3_stmt *stmt, const int position, MDObject &valueOut);

    static char *errmsg;
    static const char *zLeftover;
    static int rc;
    static sqlite3_stmt *stmt;

    static std::stringstream preparedStream;	// Stream.
    static sqlite3_stmt * preparedStmt;	// SQL statement.

    ///Non-static attributes
    int tableId;
    MetaDataDb *myMd;
    MDCache *myCache;

    friend class MDSqlStaticInit;
    friend class MetaDataDb;
    friend class MDIterator;
    ///similar to "operator"
    bool equals(const MDSql &op);

}
;//close class MDSql


/** Class to store some cached sql statements.
 */
class MDCache
{
public:
    sqlite3_stmt *iterStmt;
    std::map<MDLabel, sqlite3_stmt*> getValueCache;
    std::map<MDLabel, sqlite3_stmt*> setValueCache;
    sqlite3_stmt *addRowStmt;

    MDCache();
    ~MDCache();
    void clear();
};

/** Just to work as static constructor for initialize database.
 */
class MDSqlStaticInit
{
private:
    MDSqlStaticInit()
    {
        MDSql::sqlBegin();
    }//close constructor

    ~MDSqlStaticInit()
    {
        MDSql::sqlEnd();
    }//close destructor

    friend class MDSql;
}
;//close class MDSqlStaticInit

#endif
