/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
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

#ifndef CORE_UTILS_SQL_UTILS_H_
#define CORE_UTILS_SQL_UTILS_H_

#include <sqlite3.h>
#include "../metadata_label.h"
#include "../metadata_object.h"
#include "../metadata_static.h"

class sqlUtils {
public:
    /**
     * Add a new column to a DB table.
     * WARNING: this operation is typically 'very slow'.
     * If you can, define all columns of the table at the construction time.
     * @param columns to be added
     * @param db to be altered
     * @param table to be altered
     * @returns true on success
     */
    static bool addColumns(const std::vector<MDLabel> &columns,
            sqlite3 *db, const std::string &table);

    /**
     * Report last error registered by the db to std::cerr
     * @param db where the error occured
     */
    static void reportLastError(sqlite3 *db) {

    }

    /**
     * Retrieve all values within a specific row
     * @param rowId of the row
     * @param db to be read
     * @param table to be read
     * @param values will be stored here
     * @return true on success
    */
    static bool select(size_t rowId,
            sqlite3 *db, const std::string &table,
            std::vector<MDObject> &values);

    /**
     * Retrieve all rows from a table
     * All MDObjects are expected to be in the same order
     * @param db to be read
     * @param table to be read
     * @columns to be read
     * @param rows will be stored here
     * @return true on success
     */
    static bool select(sqlite3 *db, const std::string &table,
            const std::vector<MDObject> &columns,
            std::vector<std::vector<MDObject>> &rows);

    /**
     * Retrieve all values within a single column
     * @param label of the column
     * @param db to be read
     * @param table to be read
     * @param values will be stored here (appended to the end)
     * @return true on success
     */
    template<typename T>
    static bool select(const MDLabel &label,
            sqlite3 *db, const std::string &table,
            std::vector<T> &values) {
        // assuming all records are the same
        auto query = createSelectQuery({label}, table);

        sqlite3_stmt *stmt = nullptr;
        // FIXME currently, whole db is in one huge transaction. Finish whatever might be pending,
        // do our business in a clean transaction and start a new transaction after (not to break the original code)
        commitTrans(db);
        beginTrans(db);
        sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);

        MDObject obj(label);
        // execute, extract value from each row
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            extractValue(stmt, 0, obj);
            values.emplace_back(obj.getValue2(T()));
        }

        sqlite3_reset(stmt);
        sqlite3_finalize(stmt);
        endTrans(db);
        beginTrans(db);

        return checkError(db);
    }

    /**
     * Add multiple new rows into the table
     * @param records to be added (respective columns are expected to exists in the table)
     * @param db to be altered
     * @param table to be altered
     * @return true on success
     */
    static bool insert(const std::vector<std::vector<const MDObject*>> &records,
            sqlite3 *db, const std::string &table);

    /**
     * Add a single new row into the table
     * @param values in the row to be added (respective columns are expected to exists in the table)
     * @param db to be altered
     * @param table to be altered
     * @return true on success
     */
    static bool insert(const std::vector<const MDObject*> &values,
            sqlite3 *db, const std::string &table);

    /**
     * Update a single row in the table
     * @param values in the row to be updated (respective columns are expected to exists in the table)
     * @param db to be altered
     * @param table to be altered
     * @param id of the row
     * @return true on success
     */
    static bool update(const std::vector<const MDObject*> &values,
            sqlite3 *db, const std::string &table, size_t id);

protected:
    static inline void beginTrans(sqlite3 *db) {
        sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);
    }

    static inline void endTrans(sqlite3 *db) {
        sqlite3_exec(db, "END TRANSACTION", nullptr, nullptr, nullptr);
    }

    static inline void commitTrans(sqlite3 *db) {
        sqlite3_exec(db, "COMMIT TRANSACTION", nullptr, nullptr, nullptr);
    }

    static bool checkError(sqlite3 *db);

    /** Create a query for selecting multiple values from a table */
    static std::string createSelectQuery(
            const std::vector<MDObject> &values,
            const std::string &table);

    /** Create a query for selecting multiple values from a specific row */
    static std::string createSelectQuery(size_t rowId,
            const std::vector<MDObject> &values,
            const std::string &table);

    static std::string createInsertQuery(
            const std::vector<const MDObject*> &values,
            const std::string &table);

    /** Create an update query for setting multiple values from a specific row */
    static std::string createUpdateQuery(
        const std::vector<const MDObject*> &values,
        const std::string &table,
        size_t id);

private:
    /** FIXME this is copied directly from metadata_sql.h */
    static void extractValue(sqlite3_stmt *stmt, const int position, MDObject &valueOut);
    static int bindValue(sqlite3_stmt *stmt, const int position, const MDObject &valueIn);
};

#endif /* CORE_UTILS_SQL_UTILS_H_ */
