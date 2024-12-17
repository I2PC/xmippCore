/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include "xmipp_image_base.h"
#include "xmipp_error.h"
#include "xmipp_memory.h"

#include "metadata_static.h"
#include "multidim_array_base.h"
#include "xmipp_funcs.h"

#include <memory>
/*
        Base on rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321       Modified: 20030723
*/


#define MRCSIZE    1024 // Minimum size of the MRC header (when nsymbt = 0)

///@defgroup MRC MRC File format
///@ingroup ImageFormats

/** MRC Old Header
  * @ingroup MRC
  * see: http://www.ccpem.ac.uk/mrc_format/mrc2014.php for details
*/
struct MRCheadold
{          // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=uchar,1=short,2=float
    int nxStart;         //  4           unit cell offset
    int nyStart;         //  5
    int nzStart;         //  6
    int mx;              //  7           unit cell size in voxels
    int my;              //  8
    int mz;              //  9
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float amin;          // 19           minimum density value
    float amax;          // 20   80      maximum density value
    float amean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[25];     // 24           user-defined info
    float xOrigin;       // 49           phase origin (pixels) or origin of subvolume (A)
    float yOrigin;       // 50
    float zOrigin;       // 51
    char map[4];         // 52           character string 'MAP ' to identify file type
	int machst;          // 53           machine stamp encoding byte ordering of data
	float rms;           // 54           rms deviation of map from mean density
    int nlabl;           // 55           number of labels used
    char labels[10][80]; // 56-255       10 80-character labels
} ;

/** MRC Header
  * @ingroup MRC
*/
struct MRChead
{             // file header for MRC data
    int32_t nx;              //  0   0       image size
    int32_t ny;              //  1   4
    int32_t nz;              //  2   8
    int32_t mode;            //  3           0=char,1=short,2=float,6=uint16
    int32_t nxStart;         //  4           unit cell offset
    int32_t nyStart;         //  5
    int32_t nzStart;         //  6
    int32_t mx;              //  7           unit cell size in voxels
    int32_t my;              //  8
    int32_t mz;              //  9    1=Image or images stack, if volume mz=nz.
                         //              If ispg=401 then nz=number of volumes in stack * volume z dimension
                         //              mz=volume zdim
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int32_t mapc;            // 16           column axis
    int32_t mapr;            // 17           row axis
    int32_t maps;            // 18           section axis
    float amin;          // 19           minimum density value
    float amax;          // 20   80      maximum density value
    float amean;         // 21           average density value
    int32_t ispg;            // 22           space group number   0=Image/stack,1=Volume,401=volumes stack
    int32_t nsymbt;          // 23           bytes used for sym. ops. table
    float extra[25];     // 24           user-defined info
    float xOrigin;       // 49           phase origin in pixels FIXME: is in pixels or [L] units?
    float yOrigin;       // 50
    float zOrigin;       // 51
    char map[4];         // 52       identifier for map file ("MAP ")
    char machst[4];      // 53           machine stamp
    float arms;          // 54       RMS deviation
    int32_t nlabl;           // 55           number of labels used
    char labels[800];    // 56-255       10 80-character labels
} ;

int ImageBase::readMRC(size_t start_img, size_t batch_size, bool isStack /* = false */)
{

#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readMRC: Reading MRC file\n");
#endif

    std::unique_ptr< MRChead > header( new MRChead() );

    int errCode = 0;

    if ( fread( header.get(), MRCSIZE, 1, fimg ) < 1 )
        return(-2);

    // Determine byte order and swap bytes if from little-endian machine
    if ( (swap = (( abs( header->mode ) > SWAPTRIG ) || ( abs(header->nz) > SWAPTRIG ))) )
    {
#ifdef DEBUG
        fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
#endif

        swapPage((char *) header.get(), MRCSIZE - 800, DT_Float); // MRCSIZE - 800 is to exclude labels from swapping
    }

    // Convert VAX floating point types if necessary
    // if ( header->amin > header->amax )
    //    REPORT_ERROR(ERR_IMG_NOREAD,"readMRC: amin > max: VAX floating point conversion unsupported");

    size_t _xDim,_yDim,_zDim,_nDim;
    size_t _nDimSet;

    _xDim = header->nx;
    _yDim = header->ny;
    _zDim = header->nz;

    // Reading and storing axis order
    axisOrder[0] = 0;
    axisOrder[1] = 4 - header->maps;
    axisOrder[2] = 4 - header->mapr;
    axisOrder[3] = 4 - header->mapc;

    bool isVolStk = (header->ispg > 400);

    /* isStack is already true if file uses our customized "mrcs" extension. In this case
     * we ignore the stack behavior in header. If format is forced through ":" flag suffix,
     * then we also ignore the stack behavior in header */
    if ( !isStack && (isVolStk || !filename.contains(":")))
        isStack = ((header->ispg == 0 || isVolStk ) && (header->nsymbt == 0));

    // std::cout << "isStack = " << isStack << std::endl;

    if(isStack)
    {
        if (isVolStk)
        {
            _nDim = _zDim / header->mz;
            _zDim = header->mz;
        }
        else
        {
            _nDim = _zDim;// When isStack slices in Z are supposed to be a stack of images
            _zDim = 1;
        }

        if (batch_size == ALL_IMAGES) {
            _nDimSet = _nDim - start_img + 1;
        } else {
        _nDimSet = std::min( start_img + batch_size - 1, _nDim ) - start_img + 1;
        }


        if ( start_img > _nDim )
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readMRC: %s Image number %lu exceeds stack size %lu", this->filename.c_str(), start_img, _nDim));

    }
    else // If the reading is not like a stack, then the select_img is not taken into account and must be selected the only image
    {
        start_img = 1;
        _nDim = 1;
        _nDimSet = 1;
    }

    DataType datatype;
    switch ( header->mode )
    {
    case 0:
        datatype = DT_SChar;
        break;
    case 1:
        datatype = DT_Short;
        break;
    case 2:
        datatype = DT_Float;
        break;
    case 3:
        datatype = DT_CShort;
        break;
    case 4:
        datatype = DT_CFloat;
        break;
    case 6:
        datatype = DT_UShort;
        break;
    case 12:
        datatype = DT_HalfFloat;
        break;
    case 101:
        datatype = DT_UHalfByte;
        break;

    default:
        datatype = DT_Unknown;
        errCode = -1;
        break;
    }

    replaceNsize = _nDim;
    setDimensions(_xDim, _yDim, _zDim, _nDimSet);


    offset = MRCSIZE + header->nsymbt;
    size_t datasize_n;
    datasize_n = _xDim*_yDim*_zDim;

    // If mode is any of the fourier transforms (3,4)
    if ( header->mode > 2 && header->mode < 5 )
    {
        transform = CentHerm;
        fseek(fimg, 0, SEEK_END);
        if ( ftell(fimg) > offset + 0.8*datasize_n*gettypesize(datatype) )
            _xDim = (2 * (_xDim - 1));
        if ( header->mx%2 == 1 )
            _xDim += 1;     // Quick fix for odd x-size maps
        setDimensions(_xDim, _yDim, _zDim, _nDim);
    }

    MDMainHeader.setValue(MDL_MIN,(double)header->amin);
    MDMainHeader.setValue(MDL_MAX,(double)header->amax);
    MDMainHeader.setValue(MDL_AVG,(double)header->amean);
    MDMainHeader.setValue(MDL_STDDEV,(double)header->arms);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    double sampling;
    MDMainHeader.getValueOrDefault(MDL_SAMPLINGRATE_X, sampling, 1.0);
    if ( header->mx && header->a!=0 && sampling == 1.0)//ux
        MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double)header->a/header->mx);
    MDMainHeader.getValueOrDefault(MDL_SAMPLINGRATE_Y, sampling, 1.0);
    if ( header->my && header->b!=0 && sampling == 1.0)//yx
        MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double)header->b/header->my);
    MDMainHeader.getValueOrDefault(MDL_SAMPLINGRATE_Z, sampling, 1.0);
    if ( header->mz && header->c!=0 && sampling == 1.0)//zx
        MDMainHeader.setValue(MDL_SAMPLINGRATE_Z,(double)header->c/header->mz);

    if (dataMode==HEADER || (dataMode == _HEADER_ALL && _nDim > 1)) // Stop reading if not necessary
    {
        readData(NULL, 0, DT_Unknown, 0); // To consider axis ordering. Will not read actual data
        return errCode;
    }

    const size_t   imgStart = IMG_INDEX(start_img);
    const size_t   imgEnd = start_img + _nDimSet - 1;

    MD.clear();
    for (size_t i = 0; i < imgEnd-imgStart; i++)
        MD.push_back(std::unique_ptr<MDRowVec>(new MDRowVec(MDL::emptyHeaderVec())));

    /* As MRC does not support stacks, we use the geometry stored in the header
    for any image when we simulate the file is a stack.*/
    if (dataMode == _HEADER_ALL || dataMode == _DATA_ALL)
    {
        double aux;
        for ( size_t i = 0; i < imgEnd - imgStart; ++i )
        {
            MD[i]->setValue(MDL_SHIFT_X, (double) -header->nxStart);
            MD[i]->setValue(MDL_SHIFT_Y, (double) -header->nyStart);
            MD[i]->setValue(MDL_SHIFT_Z, (double) -header->nzStart);

            // We include auto detection of MRC2000 or CCP4 style origin based on http://situs.biomachina.org/fmap.pdf
            if (header->xOrigin != 0)
                MD[i]->setValue(MDL_ORIGIN_X, (double)-header->xOrigin);
            else if (header->nxStart != 0 && MDMainHeader.getValue(MDL_SAMPLINGRATE_X,aux))
                MD[i]->setValue(MDL_ORIGIN_X, -header->nxStart/aux);

            if (header->yOrigin !=0)
                MD[i]->setValue(MDL_ORIGIN_Y, (double)-header->yOrigin);
            else if(header->nyStart !=0 && MDMainHeader.getValue(MDL_SAMPLINGRATE_Y,aux))
                MD[i]->setValue(MDL_ORIGIN_Y, -header->nyStart/aux);

            if (header->zOrigin != 0)
                MD[i]->setValue(MDL_ORIGIN_Z, (double)-header->zOrigin);
            else if(header->nzStart !=0 && MDMainHeader.getValue(MDL_SAMPLINGRATE_Z,aux))
                MD[i]->setValue(MDL_ORIGIN_Z, -header->nzStart/aux);
        }
    }

    if ( dataMode < DATA )   // Don't read the individual header and the data if not necessary
    {
        readData(NULL, 0, DT_Unknown, 0); // To consider axis ordering. Will not read actual data
        return errCode;
    }

    // Lets read the data

    // 4-bits mode: Here is the magic to expand the compressed images
    if (datatype == DT_UHalfByte){
        readData4bit(fimg, start_img, datatype, 0);
    }
    else{
        readData(fimg, start_img, datatype, 0);
    }

    return errCode;
}

// I/O prototypes
/** MRC Reader
  * @ingroup MRC
*/
int ImageBase::readMRC(size_t select_img, bool isStack /* = false*/)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readMRC: Reading MRC file\n");
#endif

    if (select_img == ALL_IMAGES) {
        return readMRC(1, ALL_IMAGES, isStack);
    }
    return readMRC(select_img, 1, isStack);
}

/** MRC Writer
  * @ingroup MRC
*/
int ImageBase::writeMRC(size_t select_img, bool isStack, int mode, const String &bitDepth, CastWriteMode castMode)
{
    MRChead*  header = (MRChead *) askMemory(sizeof(MRChead));

    // Cast T to datatype
    DataType wDType,myTypeID = myT();

    if (bitDepth == "")
    {
        castMode = CW_CAST;
        switch(myTypeID)
        {
        case DT_Double:
        case DT_Float:
        case DT_Int:
        case DT_UInt:
            wDType = DT_Float;
            header->mode = 2;
            break;
        case DT_UShort:
            wDType = DT_UShort;
            header->mode = 6;
            break;
        case DT_Short:
            wDType = DT_Short;
            header->mode = 1;
            break;
        case DT_SChar:
            castMode = CW_CONVERT;
            /* no break */
        case DT_UChar:
            wDType = DT_UChar;
            header->mode = 0;
            break;
        case DT_CFloat:
        case DT_CDouble:
            wDType = DT_CFloat;
            header->mode = 4;
            break;
        case DT_HalfFloat:
            wDType = DT_HalfFloat;
            header->mode = 12;
            break;
        //case DT_UHalfByte:
        default:
            wDType = DT_Unknown;
            (void)wDType; // to suppress dead assignment warning
            REPORT_ERROR(ERR_TYPE_INCORRECT,(std::string)"ERROR: Unsupported data type by MRC format.");
        }
    }
    else //Convert to other data type
    {
        // Default Value
        wDType = (bitDepth == "default") ? DT_Float : datatypeRAW(bitDepth);

        switch (wDType)
        {
        case DT_Double:
        case DT_Int:
        case DT_UInt:
        case DT_Float:
            header->mode = 2;
            break;
        case DT_UChar:
            header->mode = 0;
            break;
        case DT_UShort:
            header->mode = 6;
            break;
        case DT_Short:
            header->mode = 1;
            break;
        case DT_CFloat:
        case DT_CDouble:
            header->mode = 4;
            break;
        case DT_HalfFloat:
            header->mode = 12;
            break;
        default:
            REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: incorrect MRC bits depth value.");
        }
    }

    if (mmapOnWrite)
    {
        MDMainHeader.setValue(MDL_DATATYPE,(int) wDType);
        if (!checkMmapT(wDType))
        {
            if (dataMode < DATA && castMode == CW_CAST) // This means ImageGeneric wants to know which DataType must use in mapFile2Write
                return 0;
            else //Mapping is an extra. When not available, go on and do not report an error.
            {
                /* In this case we cannot map the file because required and feasible datatypes are
                 * not compatible. Then we denote to MapFile2Write the same incoming datatype to
                 * keep using this Image object as usual, without mapping on write.
                 */
                mmapOnWrite = false;
                dataMode = DATA;
                MDMainHeader.setValue(MDL_DATATYPE,(int) myTypeID);

                // In case Image size great then, at least, map the multidimarray
                if (mdaBase->nzyxdim*gettypesize(wDType) > tiff_map_min_size)
                    mdaBase->setMmap(true);

                // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
                //if memory already allocated use it (no resize allowed)
                mdaBase->coreAllocateReuse();

                return 0;
            }
        }
        else
            dataMode = DATA;
    }


    /*
         if ( transform != NoTransform )
             img_convert_fourier(p, CentHerm);
     */

    // Map the parameters
    strncpy(header->map, "MAP ", 4);
    // FIXME TO BE DONE WITH rwCCP4!!
    //set_CCP4_machine_stamp(header->machst);
    char* machine_stamp;
    machine_stamp = (char *)(header->machst);
    if(IsLittleEndian())
    {
        machine_stamp[0] = 68;
        machine_stamp[1] = 65;
    }
    else
    {
        machine_stamp[0] = machine_stamp[1] = 17;
    }
    //                    case LittleVAX:
    //                        machine_stamp[0] = 34;
    //                        machine_stamp[1] = 65;
    //                        break;

    size_t Xdim, Ydim, Zdim, Ndim;
    getDimensions(Xdim, Ydim, Zdim, Ndim);

    /**
     * header->a,b,c info is related to sampling rate, so it is
     * only written when writing header, so it is initialized to
     * number of voxels to avoid a mistaken value.
     * If sampling is provided a, b and c are overwritten bellow
     **/

    header->mx = header->nx = Xdim;
    header->my = header->ny = Ydim;
    header->mz = header->nz = Zdim;

    // Obtaining sampling rate for each dimension and calculating cube size
    // By default sampling rate is 1.0 if no real value was found
    double sampling;
    MDMainHeader.getValueOrDefault(MDL_SAMPLINGRATE_X, sampling, 1.0);
    header->a = (float)(Xdim * sampling);
    MDMainHeader.getValueOrDefault(MDL_SAMPLINGRATE_Y, sampling, 1.0);
    header->b = (float)(Ydim * sampling);
    MDMainHeader.getValueOrDefault(MDL_SAMPLINGRATE_Z, sampling, 1.0);
    header->c = (float)(Zdim * sampling);

    if ( transform == CentHerm )
        header->nx = Xdim/2 + 1;        // If a transform, physical storage is nx/2 + 1

    header->alpha = 90.;
    header->beta  = 90.;
    header->gamma = 90.;

    //    header->mx = 0;//(int) (ua/ux + 0.5);
    //    header->my = 0;//(int) (ub/uy + 0.5);
    //    header->mz = 0;//(int) (uc/uz + 0.5);
    header->mapc = 1;
    header->mapr = 2;
    header->maps = 3;
    double aux,aux2;

    //    header->a = 0.;// ua;
    //    header->b = 0.;// ub;
    //    header->c = 0.;// uc;

    if (!MDMainHeader.empty())
    {
#define SET_MAIN_HEADER_VALUE(field, label)  MDMainHeader.getValueOrDefault(label, aux, 0.); header->field = (float)aux
        SET_MAIN_HEADER_VALUE(amin, MDL_MIN);
        SET_MAIN_HEADER_VALUE(amax, MDL_MAX);
        SET_MAIN_HEADER_VALUE(amean, MDL_AVG);
        SET_MAIN_HEADER_VALUE(arms, MDL_STDDEV);

        if ((dataMode == _HEADER_ALL || dataMode == _DATA_ALL))
        {
#define SET_HEADER_SHIFT(field, label)  MD[0]->getValueOrDefault(label, aux, 0.); header->field = -(int) round(aux)
            SET_HEADER_SHIFT(nxStart, MDL_SHIFT_X);
            SET_HEADER_SHIFT(nyStart, MDL_SHIFT_Y);
            SET_HEADER_SHIFT(nzStart, MDL_SHIFT_Z);
#define SET_HEADER_ORIGIN(field, label1, label2)  MD[0]->getValueOrDefault(label1, aux, 0.);MDMainHeader.getValueOrDefault(label2, aux2, 0.);\
              header->field = (float) (aux * aux2)

            SET_HEADER_ORIGIN(xOrigin, MDL_ORIGIN_X, MDL_SAMPLINGRATE_X);
            SET_HEADER_ORIGIN(yOrigin, MDL_ORIGIN_Y, MDL_SAMPLINGRATE_Y);
            SET_HEADER_ORIGIN(zOrigin, MDL_ORIGIN_Z, MDL_SAMPLINGRATE_Z);

#define SET_HEADER_CELL_DIM(field, label1, dimSize)  MDMainHeader.getValueOrDefault(label1, aux, 0.);\
              header->field = (float) (aux * dimSize)

            SET_HEADER_CELL_DIM(a, MDL_SAMPLINGRATE_X, Xdim);
            SET_HEADER_CELL_DIM(b, MDL_SAMPLINGRATE_Y, Ydim);
            SET_HEADER_CELL_DIM(c, MDL_SAMPLINGRATE_Z, Zdim);
        }
        else
        {
            header->nxStart = header->nyStart = header->nzStart = 0;
            header->xOrigin = header->yOrigin = header->zOrigin = 0;
        }
    }

    header->nsymbt = 0;
    header->nlabl = 10; // or zero?
    //strncpy(header->labels, p->label.c_str(), 799);

    offset = MRCSIZE + header->nsymbt;
    size_t datasize, datasize_n;
    datasize_n = Xdim*Ydim*Zdim;
    datasize = datasize_n * gettypesize(wDType);

    //#define DEBUG
#ifdef DEBUG

    printf("DEBUG rwMRC: Offset = %ld,  Datasize_n = %ld\n", offset, datasize_n);
#endif

    size_t imgStart = 0;

    if (Ndim > 1 || filename.contains(":mrcs")) // If format is forced through ":" flag suffix, then ignore the stack behavior in header
        isStack = true;

    bool isVolStk = isStack && Zdim > 1;
    size_t nDimHeader = Ndim;

    if (isStack)
    {
        imgStart = IMG_INDEX(select_img);

        if( mode == WRITE_APPEND )
        {
            imgStart = replaceNsize;
            nDimHeader = replaceNsize + Ndim;
        }
        else if( mode == WRITE_REPLACE && select_img + Ndim - 1 > replaceNsize)
        {
            nDimHeader = select_img + Ndim - 1;
        }
        //        else if (Ndim > replaceNsize)
        //            nDimHeader = Ndim;


        if (isVolStk)
        {
            header->ispg = 401;
            header->mz = Zdim;
            header->nz = Zdim * nDimHeader;
        }
        else
        {
            header->ispg = 0;
            header->nz = nDimHeader;
        }
    }
    else // To set in the header that the file is a volume not a stack
        header->ispg = (Zdim>1)? 1:0;

    //locking
    FileLock flock;
    flock.lock(fimg);

    // Write header when needed
    if(!isStack || replaceNsize < nDimHeader)
    {
        if ( swapWrite )
            swapPage((char *) header, MRCSIZE - 800, DT_Float);
        fwrite( header, MRCSIZE, 1, fimg );
    }
    freeMemory(header, sizeof(MRChead) );

    // Jump to the selected imgStart position
    fseek( fimg,offset + (datasize)*imgStart, SEEK_SET);

    size_t imgEnd = (isStack)? Ndim : 1;

    if (checkMmapT(wDType) && !mmapOnWrite && dataMode >= DATA) {
        writeData(fimg, 0, wDType, datasize_n * imgEnd, castMode);
    } else {
        for ( size_t i = 0; i < imgEnd; i++ )
        {
            // If to also write the image data or jump its size
            if (dataMode >= DATA)
            {
                if (mmapOnWrite && Ndim == 1) // Can map one image at a time only
                {
                    mappedOffset = ftell(fimg);
                    mappedSize = mappedOffset + datasize;
                    fseek(fimg, datasize-1, SEEK_CUR);
                    fputc(0, fimg);
                }
                else
                    writeData(fimg, i*datasize_n, wDType, datasize_n, castMode);
            }
            else
                fseek(fimg, datasize, SEEK_CUR);
        }    
    }

    

    // Unlock the file
    flock.unlock();

    if (mmapOnWrite)
        mmapFile();

    return(0);
}
