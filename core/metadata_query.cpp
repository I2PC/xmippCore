#include "metadata_query.h"

String MDQuery::limitString() const
{
    if (limit == -1 && offset > 0)
        REPORT_ERROR(ERR_MD_SQL, "Sqlite does not support OFFSET without LIMIT");
    std::stringstream ss;
    if (limit != -1)
        ss << " LIMIT " << limit << " ";
    if (offset > 0)
        ss << " OFFSET " << offset << " ";
    return ss.str();
}

String MDValueRange::queryStringFunc() const
{
    std::stringstream ss;
    ss << "(" << query1->queryStringFunc() << " AND " << query2->queryStringFunc() << ")";
    return ss.str();
}

String MDMultiQuery::queryStringFunc() const
{
    if (queries.size() > 0)
    {
        std::stringstream ss;
        ss << "(" << queries[0]->queryStringFunc() << ") ";
        for (size_t i = 1; i < queries.size(); i++)
            ss << operations[i] << " (" << queries[i]->queryStringFunc() << ") ";

        return ss.str();
    }
    return " ";
}
