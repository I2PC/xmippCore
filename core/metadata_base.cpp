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

MetaData::~MetaData() {}

bool MetaData::setValueFromStr(const MDLabel label, const String &value, size_t id) {
    addLabel(label);

    if (id == BAD_OBJID)
    {
        REPORT_ERROR(ERR_MD_NOACTIVE, "setValue: please provide objId other than -1");
        exit(1);
    }
    MDObject mdValue(label);
    mdValue.fromString(value);
    this->setValue(mdValue, id);
    return true;
}

bool MetaData::getStrFromValue(const MDLabel label, String &strOut, size_t id) const {
    MDObject mdValueOut(label);
    if (!getValue(mdValueOut, id))
        return false;
    strOut = mdValueOut.toString();
    return true;
}

void MetaData::removeDisabled() {
    if (containsLabel(MDL_ENABLED))
        removeObjects(MDValueLE(MDL_ENABLED, 0)); // Remove values -1 and 0 on MDL_ENABLED label
}

WriteModeMetaData metadataModeConvert(String mode) {
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

void MetaData::baseClear() {
    _path.clear();
    _comment.clear();
    _fastStringSearch.clear();
    _fastStringSearchLabel = MDL_UNDEFINED;

    _isColumnFormat = true;
    _inFile = FileName();
    _precision = 1000;
}

void MetaData::copyInfo(const MetaData& md) {
    _path = md._path;
    _comment = md._comment;
    _fastStringSearch = md._fastStringSearch;
    _fastStringSearchLabel = md._fastStringSearchLabel;

    _isColumnFormat = md._isColumnFormat;
    _inFile = md._inFile;
}

double MetaData::precision() const {
    return 1./this->_precision;
}

bool MetaData::nextBlock(mdBuffer &buffer, mdBlock &block)
{
    BLOCK_INIT(block);
    if (buffer.size == 0)
        return false;
    // Search for data_ after a newline
    block.begin = BUFFER_FIND(buffer, "data_", 5);

    if (block.begin) // data_ FOUND!!!
    {
        block.begin += 5; //Shift data_
        size_t n = block.begin - buffer.begin;
        BUFFER_MOVE(buffer, n);
        //Search for the end of line
        char *newLine = BUFFER_FIND(buffer, "\n", 1);
        //Calculate length of block name, counting after data_
        block.nameSize = newLine - buffer.begin;
        //Search for next block if exists one
        //use assign and check if not NULL at same time
        if (!(block.end = BUFFER_FIND(buffer, "\ndata_", 6))) {
            block.end = block.begin + buffer.size; }
        else {
            block.end += 1; // to include terminal \n
        }
        block.loop = BUFFER_FIND(buffer, "\nloop_", 6);
        //If loop_ is not found or is found outside block
        //scope, the block is in column format
        if (block.loop)
        {
            if (block.loop < block.end)
                block.loop += 6; // Shift \nloop_
            else
                block.loop = NULL;
        }
        //Move buffer to end of block
        n = block.end - buffer.begin;
        BUFFER_MOVE(buffer, n);
        return true;
    }

    return false;
}
