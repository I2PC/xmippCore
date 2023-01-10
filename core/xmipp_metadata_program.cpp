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

#include "xmipp_metadata_program.h"
#include "xmipp_funcs.h"
#include "xmipp_image_generic.h"
#include "metadata_extension.h"

/// Empty constructor
XmippMetadataProgram::XmippMetadataProgram()
{
    oroot = oext = fn_out = fn_in = "";
    mode = MD_OVERWRITE;
    apply_geo = false;
    allow_apply_geo = false;
    produces_an_output = false;
    produces_a_metadata = false;
    each_image_produces_an_output = false;
    allow_time_bar = true;
    decompose_stacks = true;
    delete_output_stack = true;
    get_image_info = true;
    remove_disabled = true;
    single_image = input_is_metadata = input_is_stack = output_is_stack = false;
    mdInSize = 0;
    ndimOut = zdimOut = ydimOut = xdimOut = 0;
    image_label = MDL_IMAGE;
    delete_mdIn = false;
    //Flags to store metadata when -o is stack
    save_metadata_stack = false;
    keep_input_columns = false;
    track_origin = false;
    create_empty_stackfile = false;
    datatypeOut = DT_Double;
    mdIn = nullptr;
}

XmippMetadataProgram::~XmippMetadataProgram()
{
    if (delete_mdIn)
        delete mdIn;
}

int XmippMetadataProgram::tryRead(int argc, const char ** argv, bool reportErrors)
{
    this->read( argc, argv,  reportErrors );
    return errorCode;
}

void XmippMetadataProgram::init()
{}

void XmippMetadataProgram::initComments()
{
    XmippProgram::initComments();

    CommentList comments;
    comments.addComment("Input file: metadata, stack, volume or image.");
    defaultComments["-i"]=comments;

    comments.clear();
    comments.addComment("Output file: metadata, stack, volume or image.");
    defaultComments["-o"]=comments;

    comments.clear();
    comments.addComment("Rootname of output individual images.");
    comments.addComment("Output image format can be set adding extension after rootname as \":ext\".");
    defaultComments["--oroot"]=comments;
}


void XmippMetadataProgram::defineParams()
{
    processDefaultComment("-i","-i <input_file>");
    addParamsLine("alias --input;");
    addParamsLine(" [--mode+ <mode=overwrite>]   : Metadata writing mode.");
    addParamsLine("    where <mode>");
    addParamsLine("     overwrite   : Replace the content of the file with the Metadata");
    addParamsLine("     append      : Write the Metadata as a new block, removing the old one");
    defineLabelParam();

    if (produces_a_metadata)
        produces_an_output = true;

    if (each_image_produces_an_output)
    {
        processDefaultComment("-o","[-o <output_file=\"\">]");
        addParamsLine("   alias --output;");
        processDefaultComment("--oroot","[--oroot <root=\"\">]");
    }
    else if (produces_an_output)
    {
        processDefaultComment("-o","[-o <output_file=\"\">]");
        addParamsLine("   alias --output;");
    }

    addParamsLine("  [--save_metadata_stack+ <output_md=\"\">]  : Create a metadata when the output (-o) is an stack");
    addParamsLine("                             : if --oroot is used, the metadata can be saved in -o param.");
    addParamsLine("                             : if output_md is empty, the name of the stack will be used, changing the extension to xmd.");
    addParamsLine(" [--track_origin+]   : Store the original image filename in the output ");
    addParamsLine("                     : metadata in column imageOriginal.");
    addParamsLine(" [--keep_input_columns+]   : Preserve the columns from the input metadata.");
    addParamsLine("                     : Some of the column values can be changed by the program.");

    if (allow_apply_geo)
    {
        addParamsLine("  [--dont_apply_geo]   : for 2D-images: do not apply transformation stored in metadata");
    }
}//function defineParams

void XmippMetadataProgram::defineLabelParam()
{
    addParamsLine(" [--label+ <image_label=image>]   : Label to be used to read/write images.");
}

void XmippMetadataProgram::readParams()
{
    fn_in = getParam("-i");
    mode = metadataModeConvert(getParam("--mode"));

    if (produces_an_output)
        fn_out = checkParam("-o") ? getParam("-o") : "";

    if (each_image_produces_an_output)
    {
        fn_out = checkParam("-o") ? getParam("-o") : "";
        oroot = getParam("--oroot");
    }

    if (allow_apply_geo)
        apply_geo = !checkParam("--dont_apply_geo");

    // The following flags are an "advanced" options to allow save metadata
    // when the -o is an stack, each program can define its default value
    // that's why the || construct before checkParam call
    save_metadata_stack = save_metadata_stack || checkParam("--save_metadata_stack");
    track_origin = track_origin || checkParam("--track_origin");
    keep_input_columns = keep_input_columns || checkParam("--keep_input_columns");

    MetaData * md = new MetaDataVec();
    md->read(fn_in, NULL, decompose_stacks);
    delete_mdIn = true; // Only delete mdIn when called directly from command line

    setup(md, fn_out, oroot, apply_geo, MDL::str2Label(getParam("--label")));
}//function readParams

void XmippMetadataProgram::setup(MetaData *md, const FileName &out, const FileName &oroot,
                                 bool applyGeo, MDLabel image_label)
{
    this->mdIn = md;
    this->fn_out = out;
    this->oroot = oroot;
    this->image_label = image_label;
    this->doRun = true;
    this->iter = nullptr;

    if (remove_disabled)
        mdIn->removeDisabled();

    if (mdIn->isEmpty())
        REPORT_ERROR(ERR_MD_NOOBJ, "Empty input Metadata.");

    mdInSize = mdIn->size();

    if (mdIn->isMetadataFile)
        input_is_metadata = true;
    else
    {
        if (mdInSize == 1)
            single_image = true;
        else
            input_is_stack = true;
    }

    String labelStr = MDL::label2Str(image_label);

    if (image_label == MDL_UNDEFINED)
        REPORT_ERROR(ERR_MD_BADLABEL, formatString("Unknown image label '%s'.", labelStr.c_str()));

    if (!mdIn->containsLabel(image_label))
        REPORT_ERROR(ERR_MD_MISSINGLABEL,
                     formatString("Image label '%s' is missing. See option --label.", labelStr.c_str()));

    /* Output is stack if, given a filename in fn_out, mdIn has multiple images.
     * In case no output name is given, then input is overwritten and we have to
     * check if it is stack. */
    output_is_stack = mdInSize > 1 && oroot.empty() && (!fn_out.empty() || input_is_stack);

    /* Save metadata related to output stack only if required,
     * and output is a stack.*/
    //save_metadata_stack = save_metadata_stack && output_is_stack;

    // Only delete output stack in case we are not overwriting input
    delete_output_stack = (output_is_stack && delete_output_stack) ?
                          !(fn_out.empty() && oroot.empty()) : false;

    // If the output is a stack, create empty stack file in advance to avoid concurrent access to the header
    create_empty_stackfile = (each_image_produces_an_output && output_is_stack && !fn_out.empty());

    // if create, then we need to read the dimensions of the input stack
    if (get_image_info || create_empty_stackfile)
        getImageInfo(*mdIn, xdimOut, ydimOut, zdimOut, ndimOut, datatypeOut, image_label);

    apply_geo = applyGeo;
    // if input is volume do not apply geo
    if (zdimOut > 1)
        apply_geo = false;
}//function setup

void XmippMetadataProgram::show() const
{
    if (verbose==0)
        return;
    std::cout << "Input File: " << fn_in << std::endl;
    if (apply_geo)
        std::cout << "Reading geometrical transformations stored in metadata" << std::endl;
    if (!fn_out.empty())
        std::cout << "Output File: " << fn_out << std::endl;
    if (!oroot.empty())
        std::cout << "Output Root: " << oroot << std::endl;
}

void XmippMetadataProgram::preProcess()
{}

void XmippMetadataProgram::postProcess()
{}

void XmippMetadataProgram::startProcessing()
{
    if (delete_output_stack)
        fn_out.deleteFile();

    if (create_empty_stackfile)
        createEmptyFile(fn_out, xdimOut, ydimOut, zdimOut, mdInSize, true, WRITE_OVERWRITE);

    //Show some info
    show();
    // Initialize progress bar
    time_bar_size = mdInSize;
    if (allow_time_bar && verbose && !single_image)
        init_progress_bar(time_bar_size);
    time_bar_step = CEIL((double)time_bar_size / 60.0);
    time_bar_done = 0;
}

void XmippMetadataProgram::finishProcessing()
{
    if (allow_time_bar && verbose && !single_image)
        progress_bar(time_bar_size);
    writeOutput();
}

void XmippMetadataProgram::writeOutput()
{
    if (!single_image && !mdOut.isEmpty() && !fn_out.empty())
    {
        if (produces_an_output || produces_a_metadata || !oroot.empty()) // Out as independent images
            mdOut.write(fn_out.replaceExtension("xmd"));
        else if (save_metadata_stack) // Output is stack and also save its associated metadata
        {
            FileName outFileName = getParam("--save_metadata_stack");
            if (outFileName.empty())
                outFileName = fn_out.replaceExtension("xmd");
            mdOut.write(outFileName);
        }
    }
}

void XmippMetadataProgram::showProgress()
{
    if (time_bar_step>0 && time_bar_done % time_bar_step == 0 && allow_time_bar && verbose && !single_image)
        progress_bar(time_bar_done);
}

void XmippMetadataProgram::setupRowOut(const FileName &fnImgIn, const MDRow &rowIn, const FileName &fnImgOut, MDRow &rowOut) const
{
    rowOut.clear();
    if (keep_input_columns) {
        for (const MDObject* col : rowIn)
            rowOut.setValue(*col);
    }
    rowOut.setValue(image_label, fnImgOut);
    rowOut.setValue(MDL_ENABLED, 1);

    if (track_origin)
        rowOut.setValue(MDL_IMAGE_ORIGINAL, fnImgIn);
}

void XmippMetadataProgram::wait()
{
    // In the serial implementation, we don't have to wait. This will be useful for MPI programs
}

void XmippMetadataProgram::checkPoint()
{
}

bool XmippMetadataProgram::getImageToProcess(size_t &objId, size_t &objIndex)
{
    if (nullptr == iter) {
        iter = std::unique_ptr<MetaData::id_iterator>(new MetaData::id_iterator(mdIn->ids().begin()));
        iterIndex = 0;
        time_bar_done = 0;
    } else {
        if (*iter == mdIn->ids().end()) {
            throw std::logic_error("Iterating behind the end of the metadata");
        }
        ++iterIndex;
        ++(*iter);
    }
    bool isValid = *iter != mdIn->ids().end();
    if (isValid) {
        ++time_bar_done;
        objId = **iter;
        objIndex = iterIndex;
    }
    return isValid;
}

void XmippMetadataProgram::run()
{
    FileName fnImg, fnImgOut, fullBaseName;
    mdOut.clear(); //this allows multiple runs of the same Program object

    //Perform particular preprocessing
    preProcess();

    startProcessing();

    if (!oroot.empty())
    {
        if (oext.empty())
            oext           = oroot.getFileFormat();
        oextBaseName   = oext;
        fullBaseName   = oroot.removeFileFormat();
        baseName       = fullBaseName.getBaseName();
        pathBaseName   = fullBaseName.getDir();
    }

    size_t objId;
    size_t objIndex;
    while (getImageToProcess(objId, objIndex))
    {
        ++objIndex; //increment for composing starting at 1
        auto rowIn = mdIn->getRow(objId);
        rowIn->getValue(image_label, fnImg);

        if (fnImg.empty())
            break;

        fnImgOut = fnImg;

        MDRowVec rowOut;

        if (each_image_produces_an_output)
        {
            if (!oroot.empty()) // Compose out name to save as independent images
            {
                if (oext.empty()) // If oext is still empty, then use ext of indep input images
                {
                    if (input_is_stack)
                        oextBaseName = "spi";
                    else
                        oextBaseName = fnImg.getFileFormat();
                }

                if (!baseName.empty() )
                    fnImgOut.compose(fullBaseName, objIndex, oextBaseName);
                else if (fnImg.isInStack())
                    fnImgOut.compose(pathBaseName + (fnImg.withoutExtension()).getDecomposedFileName(), objIndex, oextBaseName);
                else
                    fnImgOut = pathBaseName + fnImg.withoutExtension()+ "." + oextBaseName;
            }
            else if (!fn_out.empty() )
            {
                if (single_image)
                    fnImgOut = fn_out;
                else
                    fnImgOut.compose(objIndex, fn_out); // Compose out name to save as stacks
            }
            else
                fnImgOut = fnImg;
            setupRowOut(fnImg, *rowIn.get(), fnImgOut, rowOut);
        }
        else if (produces_a_metadata)
            setupRowOut(fnImg, *rowIn.get(), fnImgOut, rowOut);

        processImage(fnImg, fnImgOut, *rowIn.get(), rowOut);

        if (each_image_produces_an_output || produces_a_metadata)
            mdOut.addRow(rowOut);

        checkPoint();
        showProgress();
    }
    wait();

    /* Generate name to save mdOut when output are independent images. It uses as prefix
     * the dirBaseName in order not overwriting files when repeating same command on
     * different directories. If baseName is set it is used, otherwise, input name is used.
     * Then, the suffix _oext is added.*/
    if (fn_out.empty() )
    {
        if (!oroot.empty())
        {
            if (!baseName.empty() )
                fn_out = findAndReplace(pathBaseName,"/","_") + baseName + "_" + oextBaseName + ".xmd";
            else
                fn_out = findAndReplace(pathBaseName,"/","_") + fn_in.getBaseName() + "_" + oextBaseName + ".xmd";
        }
        else if (input_is_metadata) /// When nor -o neither --oroot is passed and want to overwrite input metadata
            fn_out = fn_in;
    }

    finishProcessing();

    postProcess();

    /* Reset the default values of the program in case
     * to be reused.*/
    init();
}
