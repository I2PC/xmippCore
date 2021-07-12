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

#ifndef CORE_METADATA_GENERATOR_H
#define CORE_METADATA_GENERATOR_H

#include "metadata_base.h"

/** Class to generate values for columns of a metadata*/
class MDValueGenerator {
public:
    MDLabel label; //label to which generate values

    /* Destructor*/
    virtual ~MDValueGenerator()
    {}

    /* Method to be implemented in concrete generators */
    virtual void fillValue(MetaData &md, size_t objId) = 0;
    /* Fill whole metadata */
    void fill(MetaData &md);
};//end of class MDValueGenerator

///////// Some concrete generators ////////////////
typedef enum { GTOR_UNIFORM, GTOR_GAUSSIAN, GTOR_STUDENT } RandMode;

/** MDGenerator to generate random values on columns */
class MDRandGenerator: public MDValueGenerator {
protected:
    double op1, op2, op3;
    RandMode mode;

    inline double getRandValue();
public:
    MDRandGenerator(double op1, double op2, const String &mode, double op3=0.);
    void fillValue(MetaData &md, size_t objId);
};//end of class MDRandGenerator

/** Class to fill columns with constant values */
class MDConstGenerator: public MDValueGenerator {
public:
    String value;

    MDConstGenerator(const String &value);
    void fillValue(MetaData &md, size_t objId);
};//end of class MDConstGenerator

#ifdef NEVERDEFINED
/** Class to fill columns with another metadata in row format */
class MDExpandGenerator: public MDValueGenerator {
public:
    MetaData expMd;
    FileName fn;
    MDRow row;

    void fillValue(MetaData &md, size_t objId);
};//end of class MDExpandGenerator
#endif

/** Class to fill columns with a lineal serie */
class MDLinealGenerator: public MDValueGenerator {
public:
    double initValue, step;
    size_t counter;

    MDLinealGenerator(double initial, double step);
    void fillValue(MetaData &md, size_t objId);
};//end of class MDExpandGenerator


/** @} */

#endif
