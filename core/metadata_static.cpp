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

#include <sstream>
#include "metadata_static.h"
#include "xmipp_error.h"

//This is needed for static memory allocation
MDLabelData * MDL::data[MDL_LAST_LABEL+1];
std::map<std::string, MDLabel> MDL::names;
MDLabelStaticInit MDL::initialization; //Just for initialization
MDLabel MDL::bufferIndex;

void MDL::addLabel(const MDLabel label, const MDLabelType type, const String &name, int tags)
{
  if (names.find(name) != names.end())
    REPORT_ERROR(ERR_ARG_INCORRECT, formatString("MDL::addLabel, label '%s' already exists.", name.c_str()));

    data[label] = new MDLabelData(type, name, tags);
    names[name] = label;
}

/**
 * Extra alias can be defined through the environment var XMIPP_EXTRA_ALIASES
 * the syntax is the following
 * XMIPP_EXTRA_ALIASES='anglePsi=otherAnglePsi;shiftX=otherShiftX;shiftY:otherShift'
 * The = sign will add a label, and the alias name will replace the current one (replace=True)
 * if the : sign is used, only a normal alias will added
 */
void MDL::addExtraAliases()
{
  const char * extra_aliases = getenv("XMIPP_EXTRA_ALIASES");

  if (extra_aliases)
  {
      StringVector sv, pair;
      String eq = "=", co = ":";
      tokenize(extra_aliases, sv, ";");
      MDLabel label;
      bool replace;

      for (std::vector<String>::iterator it = sv.begin(); it != sv.end(); ++it)
      {
          if (it->find(eq) != it->npos)
          {
              tokenize(*it, pair, "=");
              replace = true;
          }
          else if (it->find(co) != it->npos)
          {
              tokenize(*it, pair, co);
              replace = false;
          }
          else
              REPORT_ERROR(ERR_ARG_INCORRECT, "Invalid pair separator, use = or :");
          label = MDL::str2Label(pair[0]);
          // Add the label alias
          if (label != MDL_UNDEFINED)
              addLabelAlias(label, pair[1], replace);
          else
              REPORT_ERROR(ERR_ARG_INCORRECT,
                           formatString("Invalid label name: %s found in environment var XMIPP_EXTRA_ALIASES", pair[0].c_str()));
      }
  }
}

void MDL::addLabelAlias(const MDLabel label, const String &alias, bool replace,
                        MDLabelType type)
{
    names[alias] = label;
    if (replace)
    {
        data[label]->str = alias;
        if (type != LABEL_NOTYPE)
            data[label]->type = type;
    }
}

MDLabel MDL::getNewAlias(const String &alias, MDLabelType type)
{
    MDLabel newLabel = MDL::bufferIndex;

    if (newLabel == MDL_LAST_LABEL)
        REPORT_ERROR(ERR_ARG_INCORRECT, "Not more buffer labels to use!!!");

    addLabelAlias(newLabel, alias, true, type);
    MDL::bufferIndex = (MDLabel)((int)newLabel + 1);

    return newLabel;
}

void MDL::resetBufferIndex()
{
    MDL::bufferIndex = BUFFER_01;
}

void MDL::str2LabelVector(const String &labelsStr, std::vector<MDLabel> &labels)
{
    labels.clear();
    StringVector parts;
    splitString(labelsStr, " ", parts);
    for (size_t i = 0; i < parts.size(); ++i)
        if (MDL::isValidLabel(parts[i]))
            labels.push_back(MDL::str2Label(parts[i]));
        else
            REPORT_ERROR(ERR_PARAM_INCORRECT, formatString("Unknown label '%s' received.", parts[i].c_str()));
}

MDLabel MDL::str2Label(const String &labelName)
{
    if (names.find(labelName) == names.end())
        return MDL_UNDEFINED;
    return names[labelName];
}

String MDL::label2Str(const MDLabel &label)
{
    return (isValidLabel(label)) ? data[(int)label]->str : "";
}

String MDL::label2StrSql(const MDLabel label)
{
  String labelSqlite = "\"" + label2Str(label) + "\"";
  return labelSqlite;
}

String MDL::label2SqlColumn(const MDLabel label)
{
    std::stringstream ss;
    ss << MDL::label2StrSql(label) << " ";

    switch (MDL::labelType(label))
    {
    case LABEL_BOOL: //bools are int in sqlite3
    case LABEL_INT:
    case LABEL_SIZET:
        ss << "INTEGER";
        break;
    case LABEL_DOUBLE:
        ss << "REAL";
        break;
    case LABEL_STRING:
        ss << "TEXT";
        break;
    case LABEL_VECTOR_DOUBLE:
    case LABEL_VECTOR_SIZET:
        ss << "TEXT";
        break;
    case LABEL_NOTYPE:
        ss << "NO_TYPE";
        break;
    }
    return ss.str();
}

String MDL::labelType2Str(MDLabelType type)
{
    switch (type)
    {
    case LABEL_STRING:
        return "STRING";
    case LABEL_DOUBLE:
        return "DOUBLE";
    case LABEL_INT:
        return "INT";
    case LABEL_BOOL:
        return "BOOL";
    case LABEL_VECTOR_DOUBLE:
        return "VECTOR(DOUBLE)";
    case LABEL_SIZET:
        return "SIZE_T";
    case LABEL_VECTOR_SIZET:
        return "VECTOR(SIZE_T)";
    case LABEL_NOTYPE:
        return "NO_TYPE";
    }
    return "UNKNOWN";
}

bool MDL::isInt(const MDLabel label)
{
    return (data[label]->type == LABEL_INT);
}

bool MDL::isLong(const MDLabel label)
{
    return (data[label]->type == LABEL_SIZET);
}

bool MDL::isBool(const MDLabel label)
{
    return (data[label]->type == LABEL_BOOL);
}

bool MDL::isString(const MDLabel label)
{
    return (data[label]->type == LABEL_STRING);
}

bool MDL::isDouble(const MDLabel label)
{
    return (data[label]->type == LABEL_DOUBLE);
}

bool MDL::isVector(const MDLabel label)
{
    return (data[label]->type == LABEL_VECTOR_DOUBLE);
}

bool MDL::isVectorLong(const MDLabel label)
{
    return (data[label]->type == LABEL_VECTOR_SIZET);
}

bool MDL::isValidLabel(const MDLabel &label)
{
    return label > MDL_UNDEFINED &&
           label < MDL_LAST_LABEL &&
           data[label] != NULL;
}

bool MDL::isValidLabel(const String &labelName)
{
    return isValidLabel(str2Label(labelName));
}

MDLabelType MDL::labelType(const MDLabel label)
{
    return data[label]->type;
}

MDLabelType MDL::labelType(const String &labelName)
{
    return data[str2Label(labelName)]->type;
}

std::map<String, MDLabel>& MDL::getLabelDict()
{
    return names;
}

bool MDL::hasTag(const MDLabel label, const int tags)
{
    return data[label]->tags & tags;
}

bool MDL::isTextFile(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_TEXTFILE;
}

bool MDL::isMetadata(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_METADATA;
}

bool MDL::isCtfParam(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_CTFPARAM;
}

bool MDL::isImage(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_IMAGE;
}

bool MDL::isStack(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_STACK;
}

bool MDL::isMicrograph(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_MICROGRAPH;
}

bool MDL::isPSD(const MDLabel label)
{
    return data[label]->tags & TAGLABEL_PSD;
}

MDRowSql MDL::emptyHeaderSql() {
    MDRowSql row;
    row.resetGeo();
    row.setValue(MDL_ANGLE_ROT, 0.);
    row.setValue(MDL_ANGLE_TILT,0.);
    return row;
}

MDRowVec MDL::emptyHeaderVec() {
    MDRowVec row;
    row.resetGeo();
    row.setValue(MDL_ANGLE_ROT, 0.);
    row.setValue(MDL_ANGLE_TILT,0.);
    return row;
}

void MDL::emptifyHeader(MDRow& row) {
    row.resetGeo();
    row.setValue(MDL_ANGLE_ROT, 0.);
    row.setValue(MDL_ANGLE_TILT,0.);
}
