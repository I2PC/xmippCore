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

#include "metadata_generator.h"
#include "xmipp_funcs.h"

inline double MDRandGenerator::getRandValue() {
    switch (mode)
    {
    case GTOR_UNIFORM:
        return rnd_unif(op1, op2);
    case GTOR_GAUSSIAN:
        return rnd_gaus(op1, op2);
    case GTOR_STUDENT:
        return rnd_student_t(op3, op1, op2);
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown random type");
    }
}

MDRandGenerator::MDRandGenerator(double op1, double op2, const String &mode, double op3) {
    static bool randomized = false;

    if (!randomized)//initialize random seed just once
    {
        randomize_random_generator();
        randomized = true;
    }
    this->op1 = op1;
    this->op2 = op2;
    this->op3 = op3;
    if (mode == "uniform")
        this->mode = GTOR_UNIFORM;
    else if (mode == "gaussian")
        this->mode = GTOR_GAUSSIAN;
    else if (mode == "student")
        this->mode = GTOR_STUDENT;
    else
        REPORT_ERROR(ERR_PARAM_INCORRECT, formatString("Unknown random type '%s'", mode.c_str()));

}

void MDRandGenerator::fillValue(MetaData &md, size_t objId) {
    double aux = getRandValue();
    md.setValue(label, aux, objId);
}

MDConstGenerator::MDConstGenerator(const String &value) {
    this->value = value;
}

void MDConstGenerator::fillValue(MetaData &md, size_t objId) {
    md.setValueFromStr(label, value, objId);
}

MDLinealGenerator::MDLinealGenerator(double initial, double step) {
    this->initValue = initial;
    this->step = step;
    counter = 0;
}

void MDLinealGenerator::fillValue(MetaData &md, size_t objId) {
    double value = initValue + step * counter++;
    if (MDL::isInt(label))
        md.setValue(label, (int)value, objId);
    else if ( MDL::isLong(label))
        md.setValue(label, (size_t)value, objId);
    else
        md.setValue(label, value, objId);
}

/* Class to generate values for columns of a metadata*/
void MDValueGenerator::fill(MetaData &md) {
    for (size_t id: md.ids())
        fillValue(md, id);
}

#ifdef NEVERDEFINED
/* Class to fill columns with another metadata in row format */
void MDExpandGenerator::fillValue(MetaData &md, size_t objId) {
    if (md.getValue(label, fn, objId))
    {
        expMd.read(fn);
        if (expMd.isColumnFormat() || expMd.isEmpty())
            REPORT_ERROR(ERR_VALUE_INCORRECT, "Only can expand non empty and row formatted metadatas");
        expMd.getRow(row, expMd.firstObject());
        md.setRow(row, objId);
    }
    else
        REPORT_ERROR(ERR_MD_BADLABEL, formatString("Can't expand missing label '%s'", MDL::label2Str(label).c_str()));
}
#endif
