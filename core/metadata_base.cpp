/***************************************************************************
 *
 * Authors:     Jan Horacek (xhorace4@fi.muni.cz)
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

#include <algorithm>
#include "metadata_base.h"

MetaData::MetaData(const MetaData &md) {
    this->copyInfo(md);
}

MetaData& MetaData::operator=(const MetaData &md) {
    this->copyInfo(md);
    return *this;
}

void MetaData::copyInfo(const MetaData& md) {
    this->_fastStringSearch = md._fastStringSearch;
    this->_fastStringSearchLabel = md._fastStringSearchLabel;
    this->_path = md._path;
    this->_comment = md._comment;
    this->_isColumnFormat = md._isColumnFormat;
    this->_precision = md._precision;
    this->_inFile = md._inFile;
    this->_activeLabels = md._activeLabels;
    this->_ignoreLabels = md._ignoreLabels;
    this->_maxRows = md._maxRows;
    this->_parsedLines = md._parsedLines;
}

bool MetaData::setValueFromStr(const MDLabel label, const String &value, size_t id)
{
    addLabel(label);

    if (id == BAD_OBJID)
    {
        REPORT_ERROR(ERR_MD_NOACTIVE, "setValue: please provide objId other than -1");
        exit(1);
    }
    MDObject mdValue(label);
    mdValue.fromString(value);
    this->setValue(mdValue, id);
}

bool MetaData::getStrFromValue(const MDLabel label, String &strOut, size_t id) const
{
    MDObject mdValueOut(label);
    if (!getValue(mdValueOut, id))
        return false;
    strOut = mdValueOut.toString();
    return true;
}

WriteModeMetaData metadataModeConvert(String mode)
{
    toLower(mode);
    if (mode.npos != mode.find("overwrite"))
        return MD_OVERWRITE;
    if (mode.npos != mode.find("append"))
        return MD_APPEND;
    REPORT_ERROR(ERR_ARG_INCORRECT,"metadataModeConvert: Invalid mode");
}

bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label) {
    std::vector<MDLabel>::const_iterator location;
    location = std::find(labelsVector.begin(), labelsVector.end(), label);
    return (location != labelsVector.end());
}
