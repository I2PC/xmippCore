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

#include <fstream>
#include <algorithm>
#include "metadata_base.h"
#include "xmipp_image.h"

#include <sys/stat.h>
#include <fcntl.h>
#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif
#include "metadata_db.h"
#include "xmipp_funcs.h"

// Get the blocks available
void getBlocksInMetaDataFile(const FileName &inFile, StringVector& blockList)
{
    if (!inFile.isMetaData())
        return;
    if (inFile.getBlockName() != "")
        return;
    blockList.clear();
    String extFile = inFile.getExtension();

    if (extFile == "xml") {
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "getBlocksInMetaDataFile");
    } else if(extFile == "sqlite") {
        getBlocksInMetaDataFileDB(inFile, blockList);
    } else { //map file
        int fd;
        MetaDataDb mdAux;
        mdAux.setMaxRows(1);
        mdAux.read(inFile);
        BUFFER_CREATE(bufferMap);
        mapFile(inFile, bufferMap.begin, bufferMap.size, fd);
        BUFFER_COPY(bufferMap, buffer);
        BLOCK_CREATE(block);
        String blockName;
        while (mdAux.nextBlock(buffer, block)) {
            BLOCK_NAME(block, blockName);
            blockList.emplace_back(blockName);
        }

        unmapFile(bufferMap.begin, bufferMap.size, fd);
    }
}

// Does the blocks exist
bool existsBlockInMetaDataFile(const FileName &inFileWithBlock) {
    return existsBlockInMetaDataFile(inFileWithBlock.removeBlockName(),
                                     inFileWithBlock.getBlockName());
}

bool existsBlockInMetaDataFile(const FileName &inFile, const String& inBlock) {
    if (!inFile.isMetaData())
        return false;
    if (!inFile.getBlockName().empty())
        return inBlock == inFile.getBlockName();

    MetaDataDb MDaux(inFile);
    //map file
    int fd;
    BUFFER_CREATE(bufferMap);
    mapFile(inFile, bufferMap.begin, bufferMap.size, fd);
    BUFFER_COPY(bufferMap, buffer);
    BLOCK_CREATE(block);
    String blockName;
    bool result = false;
    while (MDaux.nextBlock(buffer, block)) {
        BLOCK_NAME(block, blockName);
        if (inBlock == blockName) {
            result = true;
            break;
        }
    }

    unmapFile(bufferMap.begin, bufferMap.size, fd);
    return result;
}

MetaData::~MetaData() {}

bool MetaData::setValueFromStr(const MDLabel label, const String &value, size_t id) {
    addLabel(label);
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
    REPORT_ERROR(ERR_ARG_INCORRECT,"metadataModeConvert: Invalid mode: "+mode);
}

bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label) {
    std::vector<MDLabel>::const_iterator location;
    location = std::find(labelsVector.begin(), labelsVector.end(), label);
    return (location != labelsVector.end());
}

void MetaData::keepLabels(const std::vector<MDLabel> &labels) {
    auto active = this->getActiveLabels();
    for (const auto &l : active) {
        if (!vectorContainsLabel(labels, l)) {
            this->removeLabel(l);
        }
    }
}

void MetaData::clear() {
    _comment.clear();
    _fastStringSearch.clear();
    _fastStringSearchLabel = MDL_UNDEFINED;

    _isColumnFormat = true;
    _inFile = FileName();
    _precision = 1000;
    _parsedLines = 0;
}

void MetaData::copyInfo(const MetaData& md) {
    _comment = md._comment;
    _fastStringSearch = md._fastStringSearch;
    _fastStringSearchLabel = md._fastStringSearchLabel;

    _isColumnFormat = md._isColumnFormat;
    _inFile = md._inFile;
}

double MetaData::precision() const {
    return 1./this->_precision;
}

bool MetaData::nextBlock(mdBuffer &buffer, mdBlock &block) {
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

/* This function will read the possible columns from the file
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels (for OLD doc files)
 */
void MetaData::_readColumns(std::istream& is, std::vector<MDObject*> & columnValues,
                            const std::vector<MDLabel>* desiredLabels) {
    String token;
    MDLabel label;

    while (is >> token)
        if (token.find('(') == String::npos)
        {
            //label is not recognized, the MDValue will be created
            //with MDL_UNDEFINED, which will be ignored while reading data
            label = MDL::str2Label(token);

            // Try to read undefined labels as String using the buffer approach
            if (label == MDL_UNDEFINED)
                label = MDL::getNewAlias(token);

            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels
            columnValues.emplace_back(new MDObject(label));
            if (label != MDL_UNDEFINED)
                addLabel(label);

        }
}

/* This function will read the possible columns from the file
 * and mark as MDL_UNDEFINED those who aren't valid labels
 * or those who appears in the IgnoreLabels vector
 * also set the activeLabels (for new STAR files)
 */
void MetaData::_readColumnsStar(mdBlock &block,
                                std::vector<MDObject*> & columnValues,
                                const std::vector<MDLabel>* desiredLabels,
                                bool addColumns,
                                size_t id) {
    char * end = block.end;
    char * newline = NULL;
    bool found_column;
    MDLabel label;
    char * iter = block.loop;
    if (!_isColumnFormat)
    {
        iter = block.begin;
        iter = END_OF_LINE() + 1; //this should point at first label, after data_XXX
    }

    do
    {
        found_column = false;
        while (iter[0] == '#') //Skip comment
            iter = END_OF_LINE() + 1;

        //trim spaces and newlines at the beginning
        while ( isspace(iter[0]))
            ++iter;

        if (iter < end && iter[0] == '_')
        {
            found_column = true;
            ++iter; //shift _
            std::stringstream ss;
            newline = END_OF_LINE();
            //Last label and no data needs this check
            if (newline == NULL)
                newline = end;
            String s(iter, newline - iter);//get current line
            ss.str(s);//set the string of the stream
            //Take the first token which is the label
            //if the label contain spaces will fail
            ss >> s; //get the first token, the label
            label = MDL::str2Label(s);
            if (label == MDL_UNDEFINED)
                label = MDL::getNewAlias(s);
            if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
                label = MDL_UNDEFINED; //ignore if not present in desiredLabels

            if (label != MDL_UNDEFINED)
                addLabel(label);

            if (addColumns)
            {
                MDObject * _mdObject = new MDObject(label);
                columnValues.emplace_back(_mdObject);//add the value here with a char
                if(!_isColumnFormat)
                    _parseObject(ss, *_mdObject, id);
            }
            iter = newline + 1;//go to next line character
        }
    }
    while (found_column)
        ;

    // This condition fails for empty blocks
    // if (iter < block.end)
    if (iter <= block.end +1)
        block.loop = iter; //Move loop pointer to position of last found column
}

/* Helper function to parse an MDObject and set its value.
 * The parsing will be from an input stream(istream)
 * and if parsing fails, an error will be raised
 */
void MetaData::_parseObject(std::istream &is, MDObject &object, size_t id) {
    object.fromStream(is);
    if (is.fail()) {
       String errorMsg = formatString("MetaData: Error parsing column '%s' value.", MDL::label2Str(object.label).c_str());
       object.failed = true;
       std::cerr << "WARNING: " << errorMsg << std::endl;
       //REPORT_ERROR(ERR_MD_BADLABEL, (String)"read: Error parsing data column, expecting " + MDL::label2Str(object.label));
    } else {
        if (object.label != MDL_UNDEFINED)
            setValue(object, id);
    }
}

/* This function parses rows data in START format
 */
void MetaData::_readRowsStar(mdBlock &block, std::vector<MDObject*> & columnValues,
                             const std::vector<MDLabel> *desiredLabels) {
    String line;
    std::stringstream ss;
    size_t n = block.end - block.loop;
    bool firstTime = true;

    if (n == 0)
        return;

    char * buffer = new char[n];
    memcpy(buffer, block.loop, n);
    char *iter = buffer, *end = iter + n, * newline = NULL;
    _parsedLines = 0; //Check how many lines the md have

    while (iter < end) { //while there are data lines
        //Adding \n position and check if NULL at the same time
        if (!(newline = END_OF_LINE()))
            newline = end;
        line.assign(iter, newline - iter);
        trim(line);

        if (!line.empty() && line[0] != '#') {
            //_maxRows would be > 0 if we only want to read some
            // rows from the md for performance reasons...
            // anyway the number of lines will be counted in _parsedLines
            if (_maxRows == 0 || _parsedLines < _maxRows) {
                std::stringstream ss(line);
                this->_parseObjects(ss, columnValues, desiredLabels, firstTime);
                firstTime = false;
            }
            _parsedLines++;
        }
        iter = newline + 1; //go to next line
    }

    delete[] buffer;
}

void MetaData::_readRows(std::istream& is, std::vector<MDObject*>& columnValues, bool useCommentAsImage) {
    String line = "";
    while (!is.eof() && !is.fail()) {
        // Move until the ';' or the first alphanumeric character
        while (is.peek() != ';' && isspace(is.peek()) && !is.eof())
            is.ignore(1);
        if (!is.eof()) {
            if (is.peek() == ';') { //is a comment
                is.ignore(1); //ignore the ';'
                getline(is, line);
                trim(line);
            } else if (!isspace(is.peek())) {
                size_t id = addObject();
                if (line != "") { //this is for old format files
                    if (!useCommentAsImage)
                        this->setValue(MDObject(MDL_STAR_COMMENT, line), id);
                    else
                        this->setValue(MDObject(MDL_IMAGE, line), id);
                }
                int nCol = columnValues.size();
                for (int i = 0; i < nCol; ++i)
                    this->_parseObject(is, *(columnValues[i]), id);
            }
        }
    }
}

void MetaData::readStar(const FileName &filename, const std::vector<MDLabel> *desiredLabels,
                        const String &blockRegExp, bool decomposeStack) {
    // First try to open the file as a metadata
    size_t id;
    FileName inFile = filename.removeBlockName();

    if (!(isMetadataFile = inFile.isMetaData())) { //if not a metadata, try to read as image or stack
        Image<char> image;
        if (decomposeStack) // If not decomposeStack it is no necessary to read the image header
            image.read(filename, HEADER);
        if (!decomposeStack || image().ndim == 1) { // single image; !decomposeStack must be first
            id = this->addObject();
            MetaData::setValue(MDL_IMAGE, filename, id);
            MetaData::setValue(MDL_ENABLED, 1, id);
        } else { // stack
            FileName fnTemp;
            for (size_t i = 1; i <= image().ndim; ++i) {
                fnTemp.compose(i, filename);
                id = addObject();
                MetaData::setValue(MDL_IMAGE, fnTemp, id);
                MetaData::setValue(MDL_ENABLED, 1, id);
            }
        }
        return;
    }

    std::ifstream is(inFile.c_str(), std::ios_base::in);
    std::stringstream ss;
    String line, token,_comment;
    std::vector<MDObject*> columnValues;

    getline(is, line); //get first line to identify the type of file

    if (is.fail())
        REPORT_ERROR(ERR_IO_NOTEXIST, formatString("MetaDataDb::read: File doesn't exists: %s", inFile.c_str()) );

    bool useCommentAsImage = false;
    this->_inFile = inFile;
    bool oldFormat = true;

    is.seekg(0, std::ios::beg);//reset the stream position to the beginning to start parsing

    if (line.find(FileNameVersion) != String::npos ||
        filename.getExtension() == "xmd" || filename.getExtension() == "star") {
        oldFormat = false;
        _comment.clear();

        // Skip comment parsing if we found the data key in the first line
        if (line.find("data_") != 0) {
            // Read comment
            //        is.ignore(256,'#');//format line
            is.ignore(256, '\n');//skip first line
            bool addspace = false;
            while (1) {
                getline(is, line);
                trim(line);
                if (line[0] == '#')
                {
                    line[0] = ' ';
                    trim(line);
                    if (addspace)
                        _comment += " " + line;
                    else
                        _comment += line;
                    addspace = true;
                } else
                    break;
            }
            this->setComment(_comment);
        }

        //map file
        int fd;
        BUFFER_CREATE(bufferMap);
        mapFile(inFile, bufferMap.begin, bufferMap.size, fd);

        BLOCK_CREATE(block);
        regex_t re;
        int rc = regcomp(&re, (blockRegExp+"$").c_str(), REG_EXTENDED|REG_NOSUB);
        if (blockRegExp.size() && rc != 0)
            REPORT_ERROR(ERR_ARG_INCORRECT, formatString("Pattern '%s' cannot be parsed: %s",
                         blockRegExp.c_str(), inFile.c_str()));
        BUFFER_COPY(bufferMap, buffer);
        bool firstBlock = true;
        bool singleBlock = blockRegExp.find_first_of(".[*+")==String::npos;

        String blockName;

        while (nextBlock(buffer, block)) {
            //startingPoint, remainingSize, firstData, secondData, firstloop))
            BLOCK_NAME(block, blockName);
            if (blockRegExp.size() == 0 || regexec(&re, blockName.c_str(), (size_t) 0, NULL, 0) == 0) {
                //Read column labels from the datablock that starts at firstData
                //Label ends at firstloop
                if ((_isColumnFormat = (block.loop != NULL))) {
                    _readColumnsStar(block, columnValues, desiredLabels, firstBlock);
                    // If block is empty, makes block.loop and block.end equal
                    if (block.loop == (block.end + 1))
                        block.loop--;
                    _readRowsStar(block, columnValues, desiredLabels);
                } else {
                    id = this->addObject();
                    _parsedLines = 1;
                    this->_readColumnsStar(block, columnValues, desiredLabels, firstBlock, id);
                }
                firstBlock = false;

                if (singleBlock)
                    break;
            }
        }

        unmapFile(bufferMap.begin, bufferMap.size, fd);
        regfree(&re);
        if (firstBlock)
            REPORT_ERROR(ERR_MD_BADBLOCK, formatString("Block: '%s': %s",
                         blockRegExp.c_str(), inFile.c_str()));

    } else if (line.find("Headerinfo columns:") != String::npos) {
        // This looks like an old DocFile, parse header
        std::cerr << "WARNING: ** You are using an old file format (DOCFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        is.ignore(256, ':'); //ignore all until ':' to start parsing column labels
        getline(is, line);
        ss.str(line);
        columnValues.emplace_back(new MDObject(MDL_UNDEFINED));
        columnValues.emplace_back(new MDObject(MDL_UNDEFINED));

        this->addLabel(MDL_IMAGE);
        this->_readColumns(ss, columnValues, desiredLabels);
        useCommentAsImage = true;

    } else {
        std::cerr << "WARNING: ** You are using an old file format (SELFILE) which is going "
        << "to be deprecated in next Xmipp release **" << std::endl;
        // I will assume that is an old SelFile, so only need to add two columns
        columnValues.emplace_back(new MDObject(MDL_IMAGE));//addLabel(MDL_IMAGE);
        columnValues.emplace_back(new MDObject(MDL_ENABLED));//addLabel(MDL_ENABLED);
    }

    if (oldFormat)
        this->_readRows(is, columnValues, useCommentAsImage);

    //free memory of column values
    int nCols = columnValues.size();
    for (int i = 0; i < nCols; ++i)
        delete columnValues[i];

    is.close();
}

void MetaData::writeStar(const FileName &outFile, const String & blockName, WriteModeMetaData mode) const {
    // Move to MetaData?
#ifdef XMIPP_MMAP
    if (outFile.hasImageExtension())
        REPORT_ERROR(ERR_IO,"MetaData:writeStar Trying to write metadata with image extension");

    struct stat file_status;
    int fd;
    char *map = NULL;
    char *tailMetadataFile = NULL; // auxiliary variable to keep metadata file tail in memory
    size_t size=-1;
    char *target, * target2 = NULL;

    //check if file exists or not block name has been given
    //in our format no two identical data_xxx strings may exists
    if (mode == MD_APPEND) {
        if (blockName.empty() || !outFile.exists())
            mode = MD_OVERWRITE;
        else {
            //does blockname exists?
            //remove it from file in this case
            // get length of file:
            if(stat(outFile.c_str(), &file_status) != 0)
                REPORT_ERROR(ERR_IO_NOPATH,"MetaData:writeStar can not get filesize for file "+outFile);
            size = file_status.st_size;
            if (size == 0)
                mode = MD_OVERWRITE;
        }

        if (mode == MD_APPEND) { //size=0 for /dev/stderr
            fd = open(outFile.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
            if (fd == -1)
                REPORT_ERROR(ERR_IO_NOPATH,"MetaData:writeStar can not read file named "+outFile);

            map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (map == MAP_FAILED)
                REPORT_ERROR(ERR_MEM_BADREQUEST,"MetaData:writeStar can not map memory ");

            // Is this a metadata formatted FILE
            if(strncmp(map,FileNameVersion.c_str(),FileNameVersion.length()) != 0) {
                mode = MD_OVERWRITE;
            } else {
                //block name
                String _szBlockName = formatString("\ndata_%s\n", blockName.c_str());
                size_t blockNameSize = _szBlockName.size();

                //search for the string
                target = (char *) _memmem(map, size, _szBlockName.c_str(), blockNameSize);

                if (target != NULL)
                {
                    target2 = (char *) _memmem(target+1, size - (target - map), "\ndata_", 6);

                    if (target2 != NULL)
                    {
                        //target block is not the last one, so we need to
                        //copy file from target2 to auxiliary memory and truncate
                        tailMetadataFile = (char *) malloc( ((map + size) - target2));
                        memmove(tailMetadataFile,target2, (map + size) - target2);
                    }
                    if (ftruncate(fd, target - map+1)==-1)  //truncate rest of the file
                        REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot truncate file");
                }
            }
            close(fd);

            if (munmap(map, size) == -1)
                REPORT_ERROR(ERR_MEM_NOTDEALLOC, "MetaData:writeStar, Can not unmap memory");
        }
    }

    std::ios_base::openmode openMode = (mode == MD_OVERWRITE) ? std::ios_base::out : std::ios_base::app;
    std::ofstream ofs(outFile.c_str(), openMode);

    write(ofs, blockName, mode);

    if (tailMetadataFile != NULL) {
        //append memory buffer to file
        //may a cat a buffer to a ofstream
        ofs.write(tailMetadataFile,(map + size) - target2);
        free(tailMetadataFile);
    }
    ofs.close();

#else

    REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif
}

void MetaData::append(const FileName &outFile) const
{
    if (outFile.exists())
    {
        std::ofstream ofs(outFile.c_str(), std::ios_base::app);
        _writeRows(ofs);
        ofs.close();
    }
    else
        write(outFile);
}
