/***************************************************************************
 *
 * Authors:    J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *             Jan Horacek (xhorace4@fi.muni.cz)
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
#include <iomanip>
#include <sstream>
#include "metadata_object.h"
#include "metadata_static.h"
#include "xmipp_error.h"
#include "xmipp_macros.h"
#include <limits>

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
    if ((type == LABEL_STRING) && (data.stringValue != nullptr))
    {
        delete data.stringValue;
        data.stringValue = nullptr;
    }
    else if ((type == LABEL_VECTOR_DOUBLE) && (data.vectorValue != nullptr))
    {
        delete data.vectorValue;
        data.vectorValue = nullptr;
    }
    else if ((type == LABEL_VECTOR_SIZET) && (data.vectorValueLong != nullptr))
    {
        delete data.vectorValueLong;
        data.vectorValueLong = nullptr;
    }

    label = obj.label;
    failed = obj.failed;
    type = obj.type;
    chr = obj.chr;

    if (type == LABEL_STRING)
        data.stringValue = new String(*(obj.data.stringValue));
    else if (type == LABEL_VECTOR_DOUBLE)
        data.vectorValue = new std::vector<double>(*(obj.data.vectorValue));
    else if (type == LABEL_VECTOR_SIZET)
        data.vectorValueLong = new std::vector<size_t>(*(obj.data.vectorValueLong));
    else
        data = obj.data;
}

MDObject::MDObject(const MDObject & obj)
{
    copy(obj);
}

MDObject & MDObject::operator = (const MDObject &obj)
{
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
MDObject::MDObject(MDLabel label, const int &v)
{
    MDOBJECT_INIT();
    this->setValue(v);
}
MDObject::MDObject(MDLabel label, const double &v)
{
    MDOBJECT_INIT();
    this->setValue(v);
}
MDObject::MDObject(MDLabel label, const bool &v)
{
    MDOBJECT_INIT();
    this->setValue(v);
}
MDObject::MDObject(MDLabel label, const String &v)
{
    MDOBJECT_INIT();
    this->data.stringValue = new String();
    this->setValue(v);
}
MDObject::MDObject(MDLabel label, const std::vector<double> &v)
{
    MDOBJECT_INIT();
    this->data.vectorValue = new std::vector<double>();
    this->setValue(v);
}
MDObject::MDObject(MDLabel label, const std::vector<size_t> &v)
{
    MDOBJECT_INIT();
    this->data.vectorValueLong = new std::vector<size_t>();
    this->setValue(v);
}
MDObject::MDObject(MDLabel label, const size_t &v)
{
    MDOBJECT_INIT();
    this->setValue(v);
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

//These getValue2 also do a compilation type checking
//when expanding templates functions and only
//will allow the supported types
//TODO: think if the type check if needed here

int& MDObject::getValue2(int) {
    labelTypeCheck(LABEL_INT);
    return this->data.intValue;
}

const int& MDObject::getValue2(int ) const {
    labelTypeCheck(LABEL_INT);
    return this->data.intValue;
}

double& MDObject::getValue2(double) {
    labelTypeCheck(LABEL_DOUBLE);
    return this->data.doubleValue;
}

const double& MDObject::getValue2(double) const {
    labelTypeCheck(LABEL_DOUBLE);
    return this->data.doubleValue;
}

bool& MDObject::getValue2(bool) {
    labelTypeCheck(LABEL_BOOL);
    return this->data.boolValue;
}

const bool& MDObject::getValue2(bool) const {
    labelTypeCheck(LABEL_BOOL);
    return this->data.boolValue;
}

String& MDObject::getValue2(String) {
    labelTypeCheck(LABEL_STRING);
    return *(this->data.stringValue);
}

const String& MDObject::getValue2(String) const {
    labelTypeCheck(LABEL_STRING);
    return *(this->data.stringValue);
}

std::vector<double>& MDObject::getValue2(std::vector<double>) {
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    return *(this->data.vectorValue);
}

const std::vector<double>& MDObject::getValue2(std::vector<double>) const {
    labelTypeCheck(LABEL_VECTOR_DOUBLE);
    return *(this->data.vectorValue);
}

std::vector<float>& MDObject::getValue2(std::vector<float>) {
    labelTypeCheck(LABEL_VECTOR_FLOAT);
    return *(this->data.vectorValueFloat);
}

const std::vector<float>& MDObject::getValue2(std::vector<float>) const {
    labelTypeCheck(LABEL_VECTOR_FLOAT);
    return *(this->data.vectorValueFloat);
}

std::vector<size_t>& MDObject::getValue2(std::vector<size_t>) {
    labelTypeCheck(LABEL_VECTOR_SIZET);
    return *(this->data.vectorValueLong);
}

const std::vector<size_t>& MDObject::getValue2(std::vector<size_t>) const {
    labelTypeCheck(LABEL_VECTOR_SIZET);
    return *(this->data.vectorValueLong);
}

size_t& MDObject::getValue2(size_t) {
    labelTypeCheck(LABEL_SIZET);
    return this->data.longintValue;
}

const size_t& MDObject::getValue2(size_t) const {
    labelTypeCheck(LABEL_SIZET);
    return this->data.longintValue;
}

float MDObject::getValue2(float) {
    return getValue2(0.); // double
}

float MDObject::getValue2(float) const {
    return getValue2(0.); // double
}

void MDObject::setValue(const int &iv)
{
    labelTypeCheck(LABEL_INT);
    this->data.intValue = iv;
}

void MDObject::setValue(const double &dv)
{
    labelTypeCheck(LABEL_DOUBLE);
    this->data.doubleValue = safeDouble(dv);
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
    const auto size = vv.size();
    this->data.vectorValue->resize(size);
    for (size_t i = 0; i < size; ++i) {
        this->data.vectorValue->operator[](i) = safeDouble(vv[i]);
    }
}

void  MDObject::setValue(const std::vector<float> &vv)
{
    labelTypeCheck(LABEL_VECTOR_FLOAT);
    const auto size = vv.size();
    this->data.vectorValueFloat->resize(size);
    for (size_t i = 0; i < size; ++i) {
        this->data.vectorValueFloat->operator[](i) = safeFloat(vv[i]);
    }
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
            {
                std::string tmp;
                is >> tmp;

                char* end = nullptr;
                data.doubleValue = std::strtod(tmp.c_str(), &end);
                
                if (end != tmp.c_str() + tmp.size())
                {
                    // Set failure flag
                    is.setstate(std::ios::failbit);
                }
            }
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
                data.vectorValue->emplace_back(d);
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
                data.vectorValueLong->emplace_back(value);
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

bool MDObject::operator==(const MDObject &obj) const {
    return this->eq(obj, 0);
}

bool MDObject::eq(const MDObject &obj, double epsilon) const {
    // FIXME: allow to compare e.g. int & double & longint
    if (this->label != obj.label)
        return false;

    if (this->type != obj.type)
        throw std::logic_error("MDObject: cannot compare == objects of different type");

    switch (this->type) {
        case LABEL_INT:
            return this->data.intValue == obj.data.intValue;

        case LABEL_BOOL:
            return this->data.boolValue == obj.data.boolValue;

        case LABEL_SIZET:
            return this->data.longintValue == obj.data.longintValue;

        case LABEL_DOUBLE:
            return std::abs(this->data.doubleValue - obj.data.doubleValue) <= epsilon;

        case LABEL_STRING:
            return *(this->data.stringValue) == *(obj.data.stringValue);

        case LABEL_VECTOR_DOUBLE:
            if (this->data.vectorValue->size() != obj.data.vectorValue->size())
                return false;
            for (size_t i = 0; i < this->data.vectorValue->size(); i++)
                if (std::abs((*this->data.vectorValue)[i] - (*obj.data.vectorValue)[i]) > epsilon)
                    return false;
            return true;

        case LABEL_VECTOR_SIZET:
            return *(this->data.vectorValueLong) == *(obj.data.vectorValueLong);

        default:
            throw std::logic_error("MDObject: unknown data type");
    };
}

bool MDObject::operator<=(const MDObject &obj) const {
    // FIXME: allow to compare e.g. int & double & longint
    if (this->type == LABEL_INT)
        return this->data.intValue <= obj.data.intValue;
    if (this->type == LABEL_BOOL)
        return this->data.boolValue <= obj.data.boolValue;
    if (this->type == LABEL_SIZET)
        return this->data.longintValue <= obj.data.longintValue;
    if (this->type == LABEL_DOUBLE)
        return this->data.doubleValue <= obj.data.doubleValue;
    if (this->type == LABEL_STRING)
        return this->data.stringValue->compare(*obj.data.stringValue) <= 0;

    throw std::logic_error("MDObject: cannot compare this type on <=");
}

bool MDObject::operator!=(const MDObject &obj) const {
    return !(*this == obj);
}

bool MDObject::operator>=(const MDObject &obj) const {
    return (!(*this <= obj) || (*this == obj));
}

bool MDObject::operator<(const MDObject &obj) const {
    return ((*this <= obj) && (*this != obj));
}

bool MDObject::operator>(const MDObject &obj) const {
    return ((*this >= obj) && (*this != obj));
}

double MDObject::safeDouble(const double v) const {
        if (std::isnan(v)) {
            // when saving NaN to sqlite3 db, the actually inserted value is sth very close to zero (0)
            // this was causing incompatibilities between (non)sqlite versions of metadata
            std::cerr << "Warning: trying to work with NaN in MDObject with label " << MDL::label2Str(this->label)
            << ". Using std::numeric_limits<double>::min() instead for backward compatilibity.\n";
            return std::numeric_limits<double>::min();
        }
        return v;
    }

float MDObject::safeFloat(const float v) const {
        if (std::isnan(v)) {
            // when saving NaN to sqlite3 db, the actually inserted value is sth very close to zero (0)
            // this was causing incompatibilities between (non)sqlite versions of metadata
            std::cerr << "Warning: trying to work with NaN in MDObject with label " << MDL::label2Str(this->label)
            << ". Using std::numeric_limits<float>::min() instead for backward compatilibity.\n";
            return std::numeric_limits<float>::min();
        }
        return v;
    }

//MDObject & MDRow::operator [](MDLabel label)
//{
//    for (iterator it = begin(); it != end(); ++it)
//        if ((*it)->label == label)
//            return *(*it);
//    MDObject * pObj = new MDObject(label);
//    emplace_back(pObj);
//
//    return *pObj;
//}
