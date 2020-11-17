/*
        Base on rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321       Modified: 20030723
*/

#ifndef CORE_RWMRC_H
#define CORE_RWMRC_H

///@defgroup MRC MRC File format
///@ingroup ImageFormats

// I/O prototypes
/** MRC Reader
  * @ingroup MRC
*/
int readMRC(size_t select_img, bool isStack = false);
int readMRC(size_t start_img, size_t batch_size, bool isStack = false);

/** MRC Writer
  * @ingroup MRC
*/
int writeMRC(size_t select_img, bool isStack=false, int mode=WRITE_OVERWRITE, const String &bitDepth="", CastWriteMode castMode = CW_CAST);

#endif
