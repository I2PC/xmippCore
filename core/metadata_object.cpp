/***************************************************************************
 *
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include <iomanip>
#include <sstream>
#include "metadata_object.h"
#include "metadata_static.h"
#include "xmipp_error.h"
#include "xmipp_macros.h"

#define DOUBLE2STREAM(d) \
        if (withFormat) {\
                (os) << std::setw(12); \
                (os) << (((d) != 0. && ABS(d) < 0.001) ? std::scientific : std::fixed);\
            } os << d;

#define INT2STREAM(i) \
        if (withFormat) os << std::setw(20); \
        os << i;
        //this must have 20 since SIZE_MAX = 18446744073709551615 size


void MDObject::copy(const MDObject &obj)
{
    label = obj.label;
    failed = obj.failed;
    type = obj.type;
    chr = obj.chr;
    if (type == LABEL_STRING)
    {
        delete data.stringValue;
        data.stringValue = new String(*(obj.data.stringValue));
    }
    else if (type == LABEL_VECTOR_DOUBLE)
    {
        delete data.vectorValue;
        data.vectorValue = new std::vector<double>(*(obj.data.vectorValue));
    }
    else if (type == LABEL_VECTOR_SIZET)
    {
        delete data.vectorValueLong;
        data.vectorValueLong = new std::vector<size_t>(*(obj.data.vectorValueLong));
    }
    else
        data = obj.data;
}

MDObject::MDObject(const MDObject & obj)
{
    data.doubleValue = 0;
    copy(obj);
}

MDObject & MDObject::operator = (const MDObject &obj)
{
    data.doubleValue = 0;
    copy(obj);
    return *this;
}

inline void MDObject::labelTypeCheck(MDLabelType checkingType) const
{
    if (this->type != checkingType)
    {
        std::stringstream ss;
        ss << "Mismatch Label (" << MDL::label2Str(label)
        << ") and value type(" << MDL::labelType2Str(checkingType) << ")";
        REPORT_ERROR(ERR_MD_BADLABEL, ss.str());
    }
}

//Just a simple constructor with the label
//don't do any type checking as have not value yet
MDObject::MDObject(MDLabel label)
{
    this->label = label;
    failed = false;
    chr = _SPACE;
    if (label != MDL_UNDEFINED)
    {
        type = MDL::labelType(label);
        if (type == LABEL_STRING)
            data.stringValue = new String;
        else if (type == LABEL_VECTOR_DOUBLE)
            data.vectorValue = new std::vector<double>;
        else if (type == LABEL_VECTOR_SIZET)
            data.vectorValueLong = new std::vector<size_t>;
    }
    else
        type = LABEL_NOTYPE;
}

/// Macro to do some basic initialization
#define MDOBJECT_INIT() this->label = label; this->type = MDL::labelType(label); this->failed = false; this->chr = _SPACE;

///Constructors for each Label supported type
///these constructor will do the labels type checking
MDObject::MDObject(MDLabel label, const int &intValue)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_INT);
    this->data.intValue = intValue;
}
MDObject::MDObject(MDLabel label, const double &doubleValue)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_DOUBLE);
    this->data.doubleValue = doubleValue;
}
MDObject::MDObject(MDLabel label, const bool &boolValue)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_BOOL);
    this->data.boolValue = boolValue;
}
MDObject::MDObject(MDLabel label, const String &stringValue)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_STRING);
    this->data.stringValue = new String(stringValue);
}
MDObject::MDObject(MDLabel label, const std::vector<double> &vectorValue)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    this->data.vectorValue = new std::vector<double>(vectorValue);
}
MDObject::MDObject(MDLabel label, const std::vector<size_t> &vectorValueLong)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_VECTOR_SIZET);
    this->data.vectorValueLong = new std::vector<size_t>(vectorValueLong);
}
MDObject::MDObject(MDLabel label, const size_t &longintValue)
{
    MDOBJECT_INIT();
    labelTypeCheck(LABEL_SIZET);
    this->data.longintValue = longintValue;

}

MDObject::~MDObject()
{
    if (type == LABEL_STRING)
        delete data.stringValue;
    else if (type == LABEL_VECTOR_DOUBLE)
        delete data.vectorValue;
    else if (type == LABEL_VECTOR_SIZET)
        delete data.vectorValueLong;
}

//These getValue also do a compilation type checking
//when expanding templates functions and only
//will allow the supported types
//TODO: think if the type check if needed here
void MDObject::getValue(int &iv) const
{
    labelTypeCheck(LABEL_INT);
    iv = this->data.intValue;
}
void MDObject::getValue(double &dv) const
{
    labelTypeCheck(LABEL_DOUBLE);
    dv = this->data.doubleValue;
}
void MDObject::getValue(bool &bv) const
{
    labelTypeCheck(LABEL_BOOL);
    bv = this->data.boolValue;
}
void MDObject::getValue(String &sv) const
{
    labelTypeCheck(LABEL_STRING);
    sv = *(this->data.stringValue);
}
void  MDObject::getValue(std::vector<double> &vv) const
{
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    vv = *(this->data.vectorValue);
}
void  MDObject::getValue(std::vector<size_t> &vv) const
{
    labelTypeCheck(LABEL_VECTOR_SIZET);
    vv = *(this->data.vectorValueLong);
}
void MDObject::getValue(size_t &lv) const
{
    labelTypeCheck(LABEL_SIZET);
    lv = this->data.longintValue;
}
void MDObject::getValue(float &floatvalue) const
{
    double tmp;
    getValue(tmp);
    floatvalue = (float) tmp;
}

void MDObject::setValue(const int &iv)
{
    labelTypeCheck(LABEL_INT);
    this->data.intValue = iv;
}
void MDObject::setValue(const double &dv)
{
    labelTypeCheck(LABEL_DOUBLE);
    this->data.doubleValue = dv;
}

void MDObject::setValue(const bool &bv)
{
    labelTypeCheck(LABEL_BOOL);
    this->data.boolValue = bv;
}

void MDObject::setValue(const String &sv)
{
    labelTypeCheck(LABEL_STRING);
    *(this->data.stringValue) = sv;
}
void  MDObject::setValue(const std::vector<double> &vv)
{
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    *(this->data.vectorValue) = vv;
}
void  MDObject::setValue(const std::vector<size_t> &vv)
{
    labelTypeCheck(LABEL_VECTOR_SIZET);
    *(this->data.vectorValueLong) = vv;
}
void MDObject::setValue(const size_t &lv)
{
    labelTypeCheck(LABEL_SIZET);
    this->data.longintValue = lv;
}
void MDObject::setValue(const float &floatvalue)
{
    setValue((double) floatvalue);
}
void MDObject::setValue(const char*  &charvalue)
{
    setValue(String(charvalue));
}

void MDObject::toStream(std::ostream &os, bool withFormat, bool isSql, bool escape) const
{
    if (label == MDL_UNDEFINED) //if undefine label, store as a literal string
        os << data.stringValue;
    else
        switch (MDL::labelType(label))
        {
        case LABEL_BOOL: //bools are int in sqlite3
            os << data.boolValue;
            break;
        case LABEL_INT:
            INT2STREAM(data.intValue);
            break;
        case LABEL_SIZET:
            INT2STREAM(data.longintValue);
            break;
        case LABEL_DOUBLE:
            DOUBLE2STREAM(data.doubleValue);
            break;
        case LABEL_STRING:
            {
                char c = _SPACE;
                if (escape)
                {
                    if (isSql || data.stringValue->find_first_of(_DQUOT) != String::npos)
                        c = _QUOT;
                    else if (data.stringValue->find_first_of(_QUOT) != String::npos)
                        c = _DQUOT;
                    else if (data.stringValue->find_first_of(_SPACE) != String::npos)
                        c = _QUOT;
                    else if (data.stringValue->empty())
                        c = _QUOT;
                }
                if (c == _SPACE)
                    os << *(data.stringValue);
                else
                    os << c << *(data.stringValue) << c;
            }
            break;
        case LABEL_VECTOR_DOUBLE:
            {
                std::vector<double> &vectorDouble = *(data.vectorValue);
                if (escape)
                    os << _QUOT << " ";
                size_t size = vectorDouble.size();
                for (size_t i = 0; i < size; i++)
                {
                    double v = vectorDouble[i];
                    DOUBLE2STREAM(v);
                    os << " ";
                }
                if (escape)
                    os << _QUOT;
            }
            break;
        case LABEL_VECTOR_SIZET:
            {
                std::vector<size_t> &vector = *(data.vectorValueLong);
                if (escape)
                    os << _QUOT << " ";
                size_t size = vector.size();
                for (size_t i = 0; i < size; i++)
                    os << vector[i] << " ";
                if (escape)
                    os << _QUOT;
            }
            break;
        case LABEL_NOTYPE:
        	if (escape) os << _QUOT;
        	os << "No type";
        	if (escape) os << _QUOT;
        	break;
        }//close switch
}//close function toStream

String MDObject::toString(bool withFormat, bool isSql) const
{
    if (type == LABEL_STRING)
    {
        return isSql ? formatString("'%s'", data.stringValue->c_str()) : *data.stringValue;
    }
    std::stringstream ss;
    toStream(ss, withFormat, isSql, isSql);

    return ss.str();
}

//bool MDValue::fromStream(std::istream &is)
std::ostream& operator<< (std::ostream& os, const MDObject &value)
{
    value.toStream(os);
    return os;
}

//bool MDValue::fromStream(std::istream &is)
std::istream& operator>> (std::istream& is, MDObject &value)
{
    value.fromStream(is);
    return is;
}

bool MDObject::fromStream(std::istream &is, bool fromString)
{
    if (label == MDL_UNDEFINED) //if undefine label, store as a literal string
    {
        String s;
        is >> s;
    }
    else
    {
        //NOTE: int, bool and long(size_t) are read as double for compatibility with old doc files
        double d;
        size_t value;
        switch (type)
        {
        case LABEL_BOOL: //bools are int in sqlite3
            is >> d;
            data.boolValue = (bool) ((int)d);
            break;
        case LABEL_INT:
            is >> d;
            data.intValue = (int) d;
            break;
        case LABEL_SIZET:
            is >> d;
            data.longintValue = (size_t) d;
            break;
        case LABEL_DOUBLE:
            is >> data.doubleValue;
            break;
        case LABEL_STRING:
            {
                data.stringValue->clear();
                String s;
                is >> s;
                char chr = s[0];
                if (chr == _QUOT || chr == _DQUOT)
                {
                    s = s.substr(1, s.size() - 1); //remove first char '
                    while (s.find_last_of(chr) == String::npos)
                    {
                        data.stringValue->append(s + " ");
                        is >> s;
                    }
                    s = s.substr(0, s.size() - 1); //remove last char '
                }
                data.stringValue->append(s);
            }
            break;
        case LABEL_VECTOR_DOUBLE:
            if (!fromString)
                is.ignore(256, _QUOT);
            //if (data.vectorValue == NULL)
            //  data.vectorValue = new std::vector<double>;
            data.vectorValue->clear();
            while (is >> d) //This will stop at ending "]"
                data.vectorValue->push_back(d);
            if (!fromString)
            {
                is.clear(); //this is for clear the fail state after found ']'
                is.ignore(256, _QUOT); //ignore the ending ']'
            }
            break;
        case LABEL_VECTOR_SIZET:
            if (!fromString)
                is.ignore(256, _QUOT);
            //if (data.vectorValue == NULL)
            //  data.vectorValue = new std::vector<double>;
            data.vectorValueLong->clear();
            while (is >> value) //This will stop at ending "]"
                data.vectorValueLong->push_back(value);
            if (!fromString)
            {
                is.clear(); //this is for clear the fail state after found ']'
                is.ignore(256, _QUOT); //ignore the ending ']'
            }
            break;
        case LABEL_NOTYPE:
        	break;
        }
    }
    return is.good();
}

bool MDObject::fromString(const String& str)
{
    if (type == LABEL_STRING)
        *data.stringValue = str;
    std::stringstream ss(str);
    return fromStream(ss, true);
}

bool MDObject::fromChar(const char * szChar)
{
    std::stringstream ss(szChar);
    return fromStream(ss);
}
//MDObject & MDRow::operator [](MDLabel label)
//{
//    for (iterator it = begin(); it != end(); ++it)
//        if ((*it)->label == label)
//            return *(*it);
//    MDObject * pObj = new MDObject(label);
//    push_back(pObj);
//
//    return *pObj;
//}