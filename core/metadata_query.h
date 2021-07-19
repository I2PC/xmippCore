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

#ifndef CORE_METADATAQUERY_H
#define CORE_METADATAQUERY_H

#include <sstream>

#include "metadata_label.h"
#include "xmipp_error.h"
#include "metadata_static.h"


/** This is the base class for queries on MetaData.
 * It is abstract, so it can not be instanciated. Queries will be very
 * helpful for performing several tasks on MetaData like importing, searching
 * or removing objects.
 */
class MDQuery
{
public:
    int limit; ///< If distint of -1 the results will be limited to this value
    int offset; ///< If distint of 0, offset elements will be discarded
    MDLabel orderLabel; ///< Label to which apply sort of the results
    bool asc;

    /** Constructor. */
    MDQuery(int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID,bool asc=true)
    {
        this->limit = limit;
        this->offset = offset;
        this->orderLabel = orderLabel;
        this->asc=asc;
    }

    /** Destructor */
    virtual ~MDQuery() {}

    /** Return the ORDER BY string to be used in SQL query */
    String orderByString() const
    {
        return (String)" ORDER BY " + MDL::label2Str(orderLabel) + (asc ? " ASC" : " DESC");
    }

    /** Return the LIMIT string to be used in SQL */
    String limitString() const;

    /** Return the WHERE string to be used in SQL query */
    String whereString() const
    {
        String queryString = this->queryStringFunc();
        return (queryString == " ") ? " " : " WHERE " + queryString + " ";
    }

    /** Return the query string, should be overrided in subclasses */
    virtual String queryStringFunc() const
    {
        return " ";
    }
}
;//End of class MDQuery

/** @} */

/** Enumeration of all posible relational queries operations */
enum RelationalOp
{
    EQ, ///< Equal
    NE, ///< Not equal
    GT, ///< Greater than
    LT, ///< Less than
    GE, ///< Greater equal
    LE ///< Less equal
};

/** Subclass of MDQuery and base for all relational queries.
 * This kind of query can be used to compare some LABEL values
 * in a column of the MetaData with an specified VALUE.
 * @see RelationalOp for possible relational operations.
 */
class MDValueRelational: public MDQuery
{
public:
    MDObject *value;
    RelationalOp op;

    template <class T>
    MDValueRelational(MDLabel label, const T &value, RelationalOp op, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        this->op = op;
        this->value = new MDObject(label, value);
    }

    MDValueRelational(const MDObject &value, RelationalOp op, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        this->op = op;
        this->value = new MDObject(value);
    }

    ~MDValueRelational()
    {
        delete this->value;
    }

    String opString() const
    {
        switch (op)
        {
        case EQ:
            return "=";
        case NE:
            return "!=";
        case GT:
            return ">";
        case LT:
            return "<";
        case GE:
            return ">=";
        case LE:
            return "<=";
        default:
        	REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown binary operator");
        }
    }

    virtual String queryStringFunc() const
    {
        return (value == NULL) ? " " : MDL::label2Str(value->label) + opString() + value->toString(false, true);
    }

    template <class T>
    void setValue(T &value)
    {
        this->value->setValue(value);
    }
}
;//end of class MDValueRelational

/** Query if MetaData column values are equal to some specific value.
 * @code
 *  ///Remove all images that are disabled
 *  MetaData md1, md2;
 *  md1.removeObjects(MDValueEQ(MDL_ENABLED, -1));
 *  ///Import objects from md2 to md1 which rot angle is 0.
 *  md1.importObjects(md2, MDValueEQ(MDL_ANGLE_ROT, 0.));
 *  @endcode
 */
class MDValueEQ: public MDValueRelational
{
public:
    template <class T>
    MDValueEQ(MDLabel label, const T &value, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, value, EQ, limit, offset, orderLabel)
    {}
}
;//end of class MDValueEQ

/** Query if MetaData column values are distint than some specific value.
 * @see MDValueEQ for examples of use.
 */
class MDValueNE: public MDValueRelational
{
public:
    template <class T>
    MDValueNE(MDLabel label, const T &value, int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, value, NE, limit, offset, orderLabel)
    {}
}
;//end of class MDValueNE

/** Query if MetaData column values are greater equal than some specific value.
 * @see MDValueEQ for examples of use.
 */
class MDValueGE: public MDValueRelational
{
public:
    template <class T>
    MDValueGE(MDLabel label, const T &valueMin, int limit = -1,int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, valueMin, GE, limit, offset, orderLabel)
    {}
}
;//end of class MDValueGE

/** Query if MetaData column values are greater than some specific value.
 * @see MDValueEQ for examples of use.
 */
class MDValueGT: public MDValueRelational
{
public:
    template <class T>
    MDValueGT(MDLabel label, const T &valueMin, int limit = -1,int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, valueMin, GT, limit, offset, orderLabel)
    {}
}
;//end of class MDValueGT

/** Query if MetaData column values are less or equal than some specific value.
 * @see MDValueEQ for examples of use.
 */
class MDValueLE: public MDValueRelational
{
public:
    template <class T>
    MDValueLE(MDLabel label, const T &valueMax, int limit = -1,int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, valueMax, LE, limit, offset, orderLabel)
    {}
}
;//end of class MDValueLE

/** Query if MetaData column values are less than some specific value.
 * @see MDValueEQ for examples of use.
 */
class MDValueLT: public MDValueRelational
{
public:
    template <class T>
    MDValueLT(MDLabel label, const T &valueMax, int limit = -1,int offset = 0, MDLabel orderLabel = MDL_OBJID)
            :MDValueRelational(label, valueMax, LT, limit, offset, orderLabel)
    {}
}
;//end of class MDValueLT

/**This subclass of Query will test if a label have a value within a minimum and maximum.
 * @code
 *  //Remove all images with rotational angle between 100 and 200
 *  MetaData md;
 *  md.removeObjects(MDValueRange(MDL_ANGLE_ROT, 100., 200.));
 *  @endcode
 */
class MDValueRange: public MDQuery
{
    MDValueRelational *query1, *query2;
public:
    MDValueRange()
    {
        query1=query2=NULL;
    }
    template <class T>
    MDValueRange(MDLabel label, const T &valueMin, const T &valueMax,
                 int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        query1 = new MDValueRelational(label, valueMin, GE);
        query2 = new MDValueRelational(label, valueMax, LE);
    }

    MDValueRange(const MDObject &o1, const MDObject &o2,
                 int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        if (o1.label != o2.label)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "Labels should be the same");
        query1 = new MDValueRelational(o1, GE);
        query2 = new MDValueRelational(o2, LE);

    }
    virtual String queryStringFunc() const;

    ~MDValueRange()
    {
        delete query1;
        delete query2;
    }
}
;//end of class MDValueRange

/**This subclass of Query will select those entries that satisfy an expression.
 * @code
 *  //Remove all images with rotational angle between 100 and 200
 *  MetaData md;
 *  md.removeObjects(MDExpression("angleRot > 100 AND angleRot < 200"));
 *  @endcode
 */
class MDExpression: public MDQuery
{
    String sExpression;
public:
    MDExpression()
    {
        sExpression = " 1=1 ";
    }
    MDExpression(String _sExpression,
                 int limit = -1,
                 int offset = 0,
                 MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        sExpression=_sExpression;
    }

    virtual String queryStringFunc() const
    {
        return sExpression;
    }

}
;//end of class MDExpression

/** Query several conditions using AND and OR.
 * This kind of query if useful if you want to check
 * two conditions at the same time, for example, import
 * all images that are enabled and have rotational angle greater than 100.
 * @code
 *  MetaData md1, md2;
 *  MDValueEQ eq(MDL_ENABLED, 1);
 *  MDValueGT gt(MDL_ANGLE_ROT, 100.);
 *  MDMultiQuery multi;
 *  //The first query added has the same effect doing with AND or OR
 *  multi.addAndQuery(eq);
 *  multi.addAndQuery(gt);
 *
 *  md1.importObjects(md2, multi);
 * @endcode
 */
class MDMultiQuery: public MDQuery
{
public:
    std::vector<const MDQuery*> queries;
    std::vector<String> operations;

    MDMultiQuery(int limit = -1, int offset = 0, MDLabel orderLabel = MDL_OBJID):MDQuery(limit, offset, orderLabel)
    {
        clear();
    }
    void addAndQuery(MDQuery &query)
    {
        queries.emplace_back(&query);
        operations.emplace_back("AND");
    }
    void addOrQuery(MDQuery &query)
    {
        queries.emplace_back(&query);
        operations.emplace_back("OR");
    }

    void clear()
    {
        queries.clear();
        operations.clear();
    }

    virtual String queryStringFunc() const;

}
;//end of class MDMultiQuery

#endif
