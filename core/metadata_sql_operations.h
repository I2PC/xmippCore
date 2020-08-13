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

#ifndef CORE_METADATA_SQL_OPERATIONS_H_
#define CORE_METADATA_SQL_OPERATIONS_H_

/** Posible Aggregation Operations in a MetaData */
enum AggregateOperation
{
    AGGR_COUNT, AGGR_MAX, AGGR_MIN, AGGR_SUM, AGGR_AVG
};

/** Posible Set Operations with MetaData */
enum SetOperation
{
    UNION, UNION_DISTINCT, INTERSECTION, SUBSTRACTION, INNER_JOIN, LEFT_JOIN, OUTER_JOIN,NATURAL_JOIN,REMOVE_DUPLICATE, DISTINCT
};

/** Enumeration of JOIN types for this operation */
enum JoinType
{
    INNER=INNER_JOIN, LEFT=LEFT_JOIN, OUTER=OUTER_JOIN,NATURAL=NATURAL_JOIN
};

#endif /* CORE_METADATA_SQL_OPERATIONS_H_ */
