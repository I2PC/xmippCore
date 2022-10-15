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

/** Union to store values */
typedef union
{
    bool boolValue;
    int intValue;
    size_t longintValue;
    double doubleValue;
    String * stringValue;
    std::vector<double> * vectorValue;
    std::vector<size_t> * vectorValueLong;
    std::vector<float> * vectorValueFloat;
} ObjectData;

/** MDObject stores single MetaData value.
 * Each column in each metadata row contains MDObject.
 * It can contain multiple types of data (see ObjectData).
 */
class MDObject
{
public:

    ObjectData data = {0};
    bool failed; // Set to True if the parsing from Star files fails
    char chr; //literal char for string, could be SPACE, QUOT or DQUOT

    void labelTypeCheck(MDLabelType checkingType) const;
    void copy(const MDObject &obj);

    MDLabel label;
    MDLabelType type = LABEL_INT;
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
    MDObject(MDLabel label, const std::vector<float> &vectorValueFloat);
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

    //These getValue2 also do a compilation type checking
    //when expanding templates functions and only
    //will allow the supported types
    //TODO: think if the type check if needed here

    // ******* WARNING ******* 
    // Methods below were orignally marked 'getValue', however they
    // took value via parameter and set value to parameter. This behavior has
    // changed. To force the programmer to read value from return type (no
    // from parameter), mehotds were temporary renamed, because compiler would
    // not otherwise fail. When whole xmipp is compiled and all occurenced of
    // old 'getValue' are replaced with 'getValue2', these methods could be
    // renamed back to 'getValue'.

    const int& getValue2(int) const;
    const double& getValue2(double) const;
    const bool&  getValue2(bool) const;
    const String& getValue2(String) const;
    const std::vector<double>& getValue2(std::vector<double>) const;
    const std::vector<float>& getValue2(std::vector<float>) const;
    const std::vector<size_t>& getValue2(std::vector<size_t>) const;
    const size_t& getValue2(size_t) const;
    float getValue2(float) const;

    int& getValue2(int);
    double& getValue2(double);
    bool&  getValue2(bool);
    String& getValue2(String);
    std::vector<double>& getValue2(std::vector<double>);
    std::vector<float>& getValue2(std::vector<float>);
    std::vector<size_t>& getValue2(std::vector<size_t>);
    size_t& getValue2(size_t);
    float getValue2(float);

    /**
     * Do not use getValue2 with char, use string.
     * chars are banned from metadata class.
     */
    void getValue2(char*) const = delete;

    void setValue(const int &iv);
    void setValue(const double &dv);
    void setValue(const bool &bv);
    void setValue(const String &sv);
    void setValue(const std::vector<double> &vv);
    void setValue(const std::vector<float> &vv);
    void setValue(const std::vector<size_t> &vv);
    void setValue(const size_t &lv);
    void setValue(const float &floatvalue);
    void setValue(const char*  &charvalue);
    void toStream(std::ostream &os, bool withFormat = false, bool isSql=false, bool escape=true) const;
    String toString(bool withFormat = false, bool isSql=false) const;
    bool fromStream(std::istream &is, bool fromString=false);
    friend std::istream& operator>> (std::istream& is, MDObject &value);
    friend std::ostream& operator<< (std::ostream& is, const MDObject &value);
    bool fromString(const String &str);
    bool fromChar(const char * str);

    bool eq(const MDObject &obj, double epsilon) const;
    bool operator == (const MDObject &obj) const;
    bool operator != (const MDObject &obj) const;
    bool operator <= (const MDObject &obj) const;
    bool operator >= (const MDObject &obj) const;
    bool operator < (const MDObject &obj) const;
    bool operator > (const MDObject &obj) const;

    friend class MDSql;

private:
    double safeDouble(const double v) const;
    float safeFloat(const float v) const;
}
; //close class MDObject

#endif
