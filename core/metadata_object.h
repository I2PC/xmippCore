/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Jan Horacek (xhorace4@fi.muni.cz)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

#ifndef CORE_METADATA_OBJECT_H
#define CORE_METADATA_OBJECT_H

#include "metadata_label.h"

#define _SPACE ' '
#define _QUOT '\''
#define _DQUOT '"'

/** Class to hold the labels values and type
 *
 */
class MDObject
{
public:

    ObjectData data;
    bool failed; // Set to True if the parsing from Star files fails
    char chr; //literal char for string, could be SPACE, QUOT or DQUOT

    void labelTypeCheck(MDLabelType checkingType) const;
    void copy(const MDObject &obj);

    MDLabel label;
    MDLabelType type;
    /** Copy constructor */
    MDObject(const MDObject & obj);
    /** Assign operator */
    MDObject & operator = (const MDObject &obj);
    //Just a simple constructor with the label
    //don't do any type checking as have not value yet
    MDObject(MDLabel label);
    ///Constructors for each Label supported type
    ///these constructor will do the labels type checking
    MDObject(MDLabel label, const int &intValue);
    MDObject(MDLabel label, const double &doubleValue);
    MDObject(MDLabel label, const bool &boolValue);
    MDObject(MDLabel label, const String &stringValue);
    MDObject(MDLabel label, const std::vector<double> &vectorValue);
    MDObject(MDLabel label, const std::vector<size_t> &vectorValueLong);
    MDObject(MDLabel label, const size_t &longintValue);

    /**
     * Do not use MDObject constructor with floats, use double.
     * Floats are banned from metadata class.
     */
    MDObject(MDLabel label, const float &floatValue) = delete;

    /**
     * Do not use MDObject constructor with char, use string.
     * Chars are banned from metadata class.
     */
    MDObject(MDLabel label, const char * &charValue) = delete;

    /// Destructor
    ~MDObject();

    //These getValue also do a compilation type checking
    //when expanding templates functions and only
    //will allow the supported types
    //TODO: think if the type check if needed here
    void  getValue(int &iv) const;
    void  getValue(double &dv) const;
    void  getValue(bool &bv) const;
    void  getValue(String &sv) const;
    void  getValue(std::vector<double> &vv) const;
    void  getValue(std::vector<size_t> &vv) const;
    void  getValue(size_t &lv) const;
    void  getValue(float &floatvalue) const;
    /**
     * Do not use getValue with char, use string.
     * chars are banned from metadata class.
     */
    void  getValue(char*  &charvalue) const = delete;

    void  setValue(const int &iv);
    void  setValue(const double &dv);
    void  setValue(const bool &bv);
    void  setValue(const String &sv);
    void  setValue(const std::vector<double> &vv);
    void  setValue(const std::vector<size_t> &vv);
    void  setValue(const size_t &lv);
    void  setValue(const float &floatvalue);
    void  setValue(const char*  &charvalue);
    void toStream(std::ostream &os, bool withFormat = false, bool isSql=false, bool escape=true) const;
    String toString(bool withFormat = false, bool isSql=false) const;
    bool fromStream(std::istream &is, bool fromString=false);
    friend std::istream& operator>> (std::istream& is, MDObject &value);
    friend std::ostream& operator<< (std::ostream& is, const MDObject &value);
    bool fromString(const String &str);
    bool fromChar(const char * str);

    friend class MDSql;
}
; //close class MDObject

#endif
