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

#include "metadata_base.h"

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
    return myMDSql->setObjectValue(id, mdValue);
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
