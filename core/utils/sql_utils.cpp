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

#include <iostream>
#include <sstream>
#include "sql_utils.h"
#include "../xmipp_error.h"

bool sqlUtils::addColumns(const std::vector<MDLabel> &columns,
        sqlite3 *db, const std::string &table) {
    // it seems that with columns, one cannot use data binding
    auto query = "ALTER TABLE " + table + " ADD COLUMN ";
    // FIXME currently, whole db is in one huge transaction. Finish whatever might be pending,
    // do our business in a clean transaction and start a new transaction after (not to break the original code)
    commitTrans(db);
    beginTrans(db);

    for (auto c : columns) {
        auto stmt = query + MDL::label2SqlColumn(c) + ";";
        sqlite3_exec(db, stmt.c_str(), nullptr, nullptr, nullptr);
    }

    endTrans(db);
    beginTrans(db);
    return checkError(db);
}

void sqlUtils::extractValue(sqlite3_stmt *stmt, const int position, MDObject &valueOut)
{
    switch (valueOut.type)
    {
    case LABEL_BOOL: //bools are int in sqlite3
        valueOut.data.boolValue = sqlite3_column_int(stmt, position) == 1;
        break;
    case LABEL_INT:
        valueOut.data.intValue = sqlite3_column_int(stmt, position);
        break;
    case LABEL_SIZET:
        valueOut.data.longintValue = sqlite3_column_int(stmt, position);
        break;
    case LABEL_DOUBLE:
        valueOut.data.doubleValue = sqlite3_column_double(stmt, position);
        break;
    case LABEL_STRING:
    {
        std::stringstream ss;
        ss << sqlite3_column_text(stmt, position);
        valueOut.data.stringValue->assign(ss.str());

        break;
    }
    case LABEL_VECTOR_DOUBLE:
    case LABEL_VECTOR_SIZET:
    {
        std::stringstream ss;
        ss << sqlite3_column_text(stmt, position);
        valueOut.fromStream(ss);
        break;
    }
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Do not know how to extract a value of type " + valueOut.type);
    }
}

bool sqlUtils::select(size_t rowId,
        sqlite3 *db, const std::string &table,
        std::vector<MDObject> &values) {
    // assuming all records are the same
    auto query = createSelectQuery(rowId, values, table);

    // FIXME currently, whole db is in one huge transaction. Finish whatever might be pending,
    // do our business in a clean transaction and start a new transaction after (not to break the original code)
    sqlite3_mutex_enter(sqlite3_db_mutex(db)); // FIXME this should be only done on demand or not at all ...
    commitTrans(db);
    beginTrans(db);

    sqlite3_stmt *stmt = nullptr;
    sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);

    // execute
    sqlite3_step(stmt);
    for (size_t i = 0; i < values.size(); ++i) {
        extractValue(stmt, i, values.at(i));
    }

    sqlite3_reset(stmt);

    sqlite3_finalize(stmt);
    endTrans(db);
    beginTrans(db);
    sqlite3_mutex_leave(sqlite3_db_mutex(db)); // FIXME this should be only done on demand or not at all ...

    return checkError(db);
}

bool sqlUtils::select(
        sqlite3 *db, const std::string &table,
        const std::vector<MDObject> &columns,
        std::vector<std::vector<MDObject>> &rows) {
    // assuming all records are the same
    auto query = createSelectQuery(columns, table);
    sqlite3_stmt *stmt = nullptr;
    // FIXME currently, whole db is in one huge transaction. Finish whatever might be pending,
    // do our business in a clean transaction and start a new transaction after (not to break the original code)
    commitTrans(db);
    beginTrans(db);
    sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);

    // execute, extract value from each row
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        rows.emplace_back(columns);
        auto &r = rows.back();
        for (size_t i = 0; i < columns.size(); ++i) {
            extractValue(stmt, i, r.at(i));
        }
    }
    sqlite3_reset(stmt);

    sqlite3_finalize(stmt);
    endTrans(db);
    beginTrans(db);

    return checkError(db);
}


std::string sqlUtils::createSelectQuery(size_t id,
        const std::vector<MDObject> &values,
        const std::string &table) {
    std::stringstream cols;
    const auto len = values.size();
    for (size_t i = 0; i < len; ++i) {
        cols << MDL::label2StrSql(values.at(i).label);
        if (len != (i + 1)) {
            cols << ", ";
        }
    }
    std::stringstream ss;
    ss << "SELECT "
        << cols.str()
        << " FROM " << table
        << " WHERE objID=" << id << ";";
    return ss.str();
}

std::string sqlUtils::createSelectQuery(
        const std::vector<MDObject> &values,
        const std::string &table) {
    std::stringstream cols;
    const auto len = values.size();
    for (size_t i = 0; i < len; ++i) {
        cols << MDL::label2StrSql(values.at(i).label);
        if (len != (i + 1)) {
            cols << ", ";
        }
    }
    std::stringstream ss;
    ss << "SELECT "
        << cols.str()
        << " FROM " << table << ";";
    return ss.str();
}

std::string sqlUtils::createUpdateQuery(
        const std::vector<const MDObject*> &values,
        const std::string &table,
        size_t id) {
    std::stringstream cols;
    const auto len = values.size();
    for (size_t i = 0; i < len; ++i) {
        cols << MDL::label2StrSql(values.at(i)->label);
        cols << "=?";
        if (len != (i + 1)) {
            cols << ", ";
        }
    }
    std::stringstream ss;
    ss << "UPDATE " << table << " SET "
            << cols.str()
            << " WHERE objID=" << id << ";";
    return ss.str();
}

bool sqlUtils::update(const std::vector<const MDObject*> &values,
            sqlite3 *db, const std::string &table, size_t id) {
    if (values.empty()) {
        return true;
    }
    // assuming all records are the same
    auto query = createUpdateQuery(values, table, id);
    sqlite3_stmt *stmt = nullptr;
    sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);
    // FIXME currently, whole db is in one huge transaction. Finish whatever might be pending,
    // do our business in a clean transaction and start a new transaction after (not to break the original code)
    commitTrans(db);
    beginTrans(db);

    // bind proper values
    for (auto i = 0; i < values.size(); ++i) {
        bindValue(stmt, i + 1, *values.at(i));
    }
    // execute
    sqlite3_step(stmt);
    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);

    sqlite3_finalize(stmt);
    endTrans(db);
    beginTrans(db);

    return checkError(db);
}

std::string sqlUtils::createInsertQuery(
        const std::vector<const MDObject*> &values,
        const std::string &table) {
    std::stringstream cols;
    std::stringstream vals;
    const auto len = values.size();
    for (size_t i = 0; i < len; ++i) {
        cols << MDL::label2StrSql(values.at(i)->label);
        vals << "?";
        if (len != (i + 1)) {
            cols << ", ";
           vals << ", ";
        }
    }
    std::stringstream ss;
    ss << "INSERT INTO " << table
            << " (" << cols.str() << ")"
            << " VALUES "
            << " (" << vals.str() << ");";
    return ss.str();
}

bool sqlUtils::insert(const std::vector<std::vector<const MDObject*>> &records,
        sqlite3 *db, const std::string &table) {
    if (0 == records.size()) {
        return true;
    }
    // assuming all records are the same
    const auto &rec = records.at(0);
    auto query = createInsertQuery(rec, table);
    sqlite3_stmt *stmt = nullptr;
    sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);
    // FIXME currently, whole db is in one huge transaction. Finish whatever might be pending,
    // do our business in a clean transaction and start a new transaction after (not to break the original code)
    commitTrans(db);
    beginTrans(db);

    const auto len = rec.size();
    for (const auto &r : records) {
        // bind proper values
        for (size_t i = 0; i < len; ++i) {
            bindValue(stmt, i + 1, *r.at(i));
        }
        // execute
        sqlite3_step(stmt);
        sqlite3_clear_bindings(stmt);
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);
    endTrans(db);
    beginTrans(db);

    return checkError(db);
}

bool sqlUtils::insert(const std::vector<const MDObject*> &values,
        sqlite3 *db, const std::string &table) {
    return insert(std::vector<std::vector<const MDObject*>>{values}, db, table);
}

int sqlUtils::bindValue(sqlite3_stmt *stmt, const int position, const MDObject &valueIn)
{
    //First reset the statement
    //rc  = sqlite3_reset(stmt);
    //std::cerr << "rc after reset: " << rc <<std::endl;
  if (valueIn.failed)
  {
    // If a value was wronly parsed, set NULL in their sqlite entry
    std::cerr << "WARNING!!! valueIn.failed = True, binding NULL" << std::endl;

    return sqlite3_bind_null(stmt, position);
  }
  else
  {
    switch (valueIn.type)
    {
    case LABEL_BOOL: //bools are int in sqlite3
        return sqlite3_bind_int(stmt, position, valueIn.data.boolValue ? 1 : 0);
    case LABEL_INT:
        return sqlite3_bind_int(stmt, position, valueIn.data.intValue);
    case LABEL_SIZET:
        return sqlite3_bind_int(stmt, position, valueIn.data.longintValue);
    case LABEL_DOUBLE:
        return sqlite3_bind_double(stmt, position, valueIn.data.doubleValue);
    case LABEL_STRING:
        return sqlite3_bind_text(stmt, position, valueIn.data.stringValue->c_str(), -1, SQLITE_TRANSIENT);
    case LABEL_VECTOR_DOUBLE:
    case LABEL_VECTOR_SIZET:
        return sqlite3_bind_text(stmt, position, valueIn.toString(false, true).c_str(), -1, SQLITE_TRANSIENT);
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Do not know how to handle this type");
    }
  }
}

bool sqlUtils::checkError(sqlite3 *db) {
        auto err = sqlite3_errcode(db);
        if (0 != err 
            && 100 != err && 101 != err) { // probably not an error in case of the sqlite3_step (fingers crossed)
            auto msg = sqlite3_errmsg(db);
            std::cerr << "SQLite3 error: " << err
                << "\n"
                << msg << std::endl;
            return false;
        }

        return true;
    }
