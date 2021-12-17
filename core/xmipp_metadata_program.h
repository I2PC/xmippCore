/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
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

#ifndef CORE_XMIPP_METADATA_PROGRAM_H_
#define CORE_XMIPP_METADATA_PROGRAM_H_

#include "xmipp_program.h"
#include "xmipp_datatype.h"
#include "xmipp_filename.h"
#include "metadata_label.h"
#include "metadata_row_sql.h"
#include "metadata_writemode.h"
#include "metadata_base.h"
#include "metadata_vec.h"


/** Special class of XmippProgram that performs some operation related with processing images.
 * It can receive a file with images(MetaData) or a single image.
 * The function processImage is virtual here and needs to be implemented by derived classes.
 * Optionally can be implemented preProcess and postProcess to perform some customs actions.
 */
class XmippMetadataProgram: public virtual XmippProgram
{
private:
    /// Input and output metadatas
    MetaData * mdIn = nullptr;
    MetaDataVec mdOut; //TODO: can be treated by reference as mdIn for
    // uses from another programs...
    std::unique_ptr<MetaData::id_iterator> iter;
    size_t iterIndex;

public:
    /// The input metadata should not be used
    /// if there is a very very special case
    /// you can use this function
    MetaData * getInputMd() { return mdIn; }
    MetaDataVec& getOutputMd() { return mdOut; }

public:
    //Image<double>   img;
    /// Filenames of input and output Metadata
    FileName fn_in, fn_out, baseName, pathBaseName, oextBaseName;
    /// Apply geo
    bool apply_geo;
    /// Output dimensions
    size_t ndimOut, zdimOut, ydimOut, xdimOut;
    DataType datatypeOut;
    /// Number of input elements
    size_t mdInSize;

protected:
    /// Metadata writing mode: OVERWRITE, APPEND
    WriteModeMetaData mode;

    /// Output extension and root
    FileName oext, oroot;
    /// MDLabel to be used to read/write images, usually will be MDL_IMAGE
    MDLabel image_label;


    // BEHAVIOR CONTROL FLAGS //

    /// Indicate that a unique final output is produced
    bool produces_an_output; // Default false (only -o param is used)
    /// Indicate that the unique final output file is a Metadata
    bool produces_a_metadata; // Default false (if true, then produces_an_output is set true)
    /// Indicate that an output is produced for each image in the input
    bool each_image_produces_an_output; // Default false (both -o --oroot params are used)
    /// Provide the program with the param --dont_apply_geo to allow the user deciding whether
    /// or not applying the transformation info as stored in the input metadata
    bool allow_apply_geo; // Default false
    /// Input Metadata will treat a stack file as a set of images instead of a unique file
    bool decompose_stacks; // Default true
    /// Delete previous output stack file prior to process images
    bool delete_output_stack; // Default true
    /// Get the input image file  dimensions to further operations
    bool get_image_info; // Default true
    /// Save the associated output metadata when output file is a stack
    bool save_metadata_stack; // Default false
    /// Include the original input image filename in the output stack
    bool track_origin; // Default false
    /// Keep input metadata columns
    bool keep_input_columns; // Default false
    /// Remove disabled images from the input selfile
    bool remove_disabled; // Default true
    /// Show process time bar
    bool allow_time_bar; // Default true

    // DEDUCED FLAGS
    /// Input is a metadata
    bool input_is_metadata;
    /// Input is a single image
    bool single_image;
    /// Input is a stack
    bool input_is_stack;
    /// Output is a stack
    bool output_is_stack;
    // Create empty output stack file prior to process images
    bool create_empty_stackfile; //
    //check whether to delete or not the input metadata
    bool delete_mdIn;

    /// Some time bar related counters
    size_t time_bar_step, time_bar_size, time_bar_done;

    virtual void initComments();
    virtual void defineParams();
    virtual void readParams();
    virtual void preProcess();
    virtual void postProcess();
    virtual void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut) = 0;
    virtual bool getImageToProcess(size_t &objId, size_t &objIndex);
    void show() const override;
    /** Do some stuff before starting processing
     * in a parallel environment usually this only be executed
     * by master.
     */
    virtual void startProcessing();
    virtual void finishProcessing();
    virtual void writeOutput(); // maybe used by checkpoint
    virtual void showProgress();

    /** Define the label param */
    virtual void defineLabelParam();

public:
    XmippMetadataProgram();

    /** Call the read function inside a try/catch block
     * The function will return the error code when
     * 0 means success
     * and a value greater than 0 represents the error type
     * */
    virtual int tryRead(int argc, const char ** argv, bool reportErrors = true);

    /** Initialization of variables should be done here
     */
    virtual void init();

    /** Setup of common XmippMetadataProgram arguments
     *  to be called from another program.
     */
    virtual void setup(MetaData *md, const FileName &o="", const FileName &oroot="",
                       bool applyGeo=false, MDLabel label=MDL_IMAGE);


    /** Destructor
     */
    virtual ~XmippMetadataProgram();

    void setMode(WriteModeMetaData _mode)
    {
        mode = _mode;
    }

    /// Prepare rowout
    void setupRowOut(const FileName &fnImgIn, const MDRow &rowIn, const FileName &fnImgOut, MDRow &rowOut) const;

    /// Wait for the distributor to finish
    virtual void wait();

    /// For very long programs, it may be needed to write checkpoints
    virtual void checkPoint();

    /// Run over all images
    virtual void run();
}
;// end of class XmippMetadataProgram

#endif /* CORE_XMIPP_METADATA_PROGRAM_H_ */
