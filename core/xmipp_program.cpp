/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (josem@cnb.csic.es)
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

#include "xmipp_program.h"
#include "argsparser.h"
#include "xmipp_program_sql.h"
#include "argsprinter.h"
#include "xmipp_funcs.h"
#include "args.h"

void XmippProgram::initComments()
{
    CommentList comments;
    comments.addComment("Verbosity level, 0 means no output.");
    defaultComments["-v"] = comments;
}

void XmippProgram::processDefaultComment(const char *param, const char *left)
{
    addParamsLine(((String)left+":"+defaultComments[param].comments[0]).c_str());
    int imax=defaultComments[param].comments.size();
    for (int i=1; i<imax; ++i)
        addParamsLine(((String)":"+defaultComments[param].comments[i]).c_str());
}

void XmippProgram::setDefaultComment(const char *param, const char *comment)
{
    defaultComments[param].clear();
    defaultComments[param].addComment(comment);
}

void XmippProgram::defineCommons()
{
    ///Add some common definitions to all Xmipp programs
    addParamsLine("== Common options ==");
    processDefaultComment("-v","[-v+ <verbose_level=1>]");
    addParamsLine("alias --verbose;");
    addParamsLine("[-h+* <param=\"\">]      : If not param is supplied show this help message.");
    addParamsLine("                         : Otherwise, specific param help is showed,");
    addParamsLine("                         : param should be provided without the '-'");
    addParamsLine("alias --help;");
    addParamsLine("[--more*]                : Show additional options.");

    ///This are a set of internal command for MetaProgram usage
    ///they should be hidden
    addParamsLine("==+++++ Internal section ==");
    addParamsLine("[--xmipp_write_definition* <dbname>] : Print metadata info about the program to sqlite database");
    addParamsLine("[--xmipp_write_wiki* ] : Print metadata info about the program in wiki format");
    addParamsLine("[--xmipp_write_autocomplete* <scriptfile>] : Add program autocomplete bash options to script file");
    addParamsLine("[--xmipp_protocol_script <script>] : This is only meanful when execute throught protocols");
    addParamsLine("[--xmipp_validate_params] : Validate input params");
}

void XmippProgram::init()
{
    initComments();
    progDef = new ProgramDef();
    this->defineParams();
    this->defineCommons();
    progDef->parse();
}

bool XmippProgram::checkBuiltIns()
{
    ///If -more_options provided, show extended usage
    if (checkParam("--more"))
        usage(1);
    ///If help requested, print usage message
    else if (checkParam("--help"))
    {
        String helpParam = getParam("-h");
        if (helpParam != "")
        {
            String cmdHelp("-");
            cmdHelp += helpParam;
            if (existsParam(cmdHelp.c_str()))
                usage(cmdHelp);
            else
            {
                cmdHelp.insert(0, "-");
                if (existsParam(cmdHelp.c_str()))
                    usage(cmdHelp);
                else
                {
                    if (verbose)
                        std::cerr << "Unrecognized param " << helpParam << " neither - or --" << std::endl;
                    usage();
                }
            }
        }
        else
            usage();
    }
    else if (checkParam("--xmipp_write_definition"))
        writeToDB();
    else if (checkParam("--xmipp_write_wiki"))
        createWiki();
    else if (checkParam("--xmipp_write_autocomplete"))
        writeToAutocomplete();
    else
        return false;
    return true;
}

void XmippProgram::writeToDB()
{
    ProgramDb db;
    db.printProgram(*progDef);
}

void XmippProgram::writeToAutocomplete( )
{
    String scriptfile = getParam("--xmipp_write_autocomplete");
    AutocompletePrinter ap(scriptfile.c_str());
    ap.printProgram(*progDef, 3);
}

void XmippProgram::createWiki()
{
    WikiPrinter wiki;
    wiki.printProgram(*progDef, 3);
}

XmippProgram::XmippProgram()
{
    //by defaul all programs have verbose = 1
    // this can be changed on mpi slaves node for no output at all
    verbose = 1;
    progDef = NULL;
    runWithoutArgs = doRun = false;
    errorCode = 0;
}

XmippProgram::XmippProgram(int argc, const char ** argv)
{
    runWithoutArgs = doRun = false;
    errorCode = 0;
    init();
    read(argc, argv);
}

XmippProgram::~XmippProgram()
{
    delete progDef;
}

void XmippProgram::defineParams()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'defineParams'");
}

void XmippProgram::run()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'run'");
}

void XmippProgram::quit(int exit_code) const
{
    exit(exit_code);
}

void XmippProgram::readParams()
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "function 'readParams'");
}



void XmippProgram::read(int argc, const char ** argv, bool reportErrors)
{
    if (progDef == NULL)
        init();

    setProgramName(argv[0]);

    doRun = false;
    errorCode = 0; //suppose no errors
    ///If not arguments are provided show the console program help
    //this behavior will be defined with environment variable XMIPP_BEHAVIOR
    if (argc == 1)
    {
        if (runWithoutArgs) {
            doRun = true;
        } else {
            usage();
        }
    }
    else
    {
       
        this->argc = argc;
        this->argv = argv;
        progDef->read(argc, argv, reportErrors);
        if (!checkBuiltIns())
        {
            if (verbose) //if 0, ignore the parameter, useful for mpi programs
                verbose = getIntParam("--verbose");
            this->readParams();
            doRun = !checkParam("--xmipp_validate_params"); //just validation, not run
        }
    }
}

void XmippProgram::read(int argc, char ** argv, bool reportErrors)
{
    read(argc,(const char **)argv,reportErrors);
}

void XmippProgram::read(const String &argumentsLine)
{
    int argc;
    char ** argv=NULL;
    char * copy=NULL;

    generateCommandLine(argumentsLine, argc, argv, copy);
    read(argc, (const char **)argv);
    delete[] copy;
    delete[] argv[0]; // the only one allocated, the rest is pointing to 'copy'
    delete[] argv;
}

int XmippProgram::tryRun()
{
    if (doRun)
        this->run();
    return errorCode;
}
/** Init progress */
void XmippProgram::initProgress(size_t total, size_t stepBin)
{
    if (verbose)
    {
        progressTotal = total;
        progressStep = XMIPP_MAX(1, total / stepBin);
        progressLast = 0;
        init_progress_bar(total);
    }
}

/** Notify progress on work */
void XmippProgram::setProgress(size_t value)
{
    progressLast = value ? value : progressLast + 1;
    if (verbose && progressLast % progressStep == 0)
        progress_bar(progressLast);
}

/** Notify end of work */
void XmippProgram::endProgress()
{
    if (verbose)
        progress_bar(progressTotal);
}

void XmippProgram::setProgramName(const char * name)
{
    progDef->name = name;
}

void XmippProgram::addUsageLine(const char * line, bool verbatim)
{
    progDef->usageComments.addComment(line,verbatim);
}
void XmippProgram::addExampleLine(const char * example, bool verbatim)
{
    progDef->examples.addComment(example, verbatim);
}
void XmippProgram::addSeeAlsoLine(const char * seeAlso)
{
    if (progDef->seeAlso=="")
        progDef->seeAlso = seeAlso;
    else
    {
        progDef->seeAlso +=", ";
        progDef->seeAlso +=seeAlso;
    }
}

void XmippProgram::clearUsage()
{
    progDef->usageComments.clear();
}
void XmippProgram::addParamsLine(const String &line)
{
    progDef->pLexer->addLine(line);
}

void XmippProgram::addParamsLine(const char * line)
{
    progDef->pLexer->addLine((String)line);
}

void XmippProgram::addKeywords(const char * keywords)
{
    progDef->keywords += " ";
    progDef->keywords += keywords;
}

const char * XmippProgram::getParam(const char * param, int arg)
{
    return progDef->getParam(param, arg);
}

const char * XmippProgram::getParam(const char * param, const char * subparam, int arg)
{
    return progDef->getParam(param, subparam, arg);
}

int XmippProgram::getIntParam(const char * param, int arg)
{
    return textToInteger(progDef->getParam(param, arg));
}

int XmippProgram::getIntParam(const char * param, const char * subparam, int arg)
{
    return textToInteger(progDef->getParam(param, subparam, arg));
}

double XmippProgram::getDoubleParam(const char * param, int arg)
{
    return textToFloat(progDef->getParam(param, arg));
}

double XmippProgram::getDoubleParam(const char * param, const char * subparam, int arg)
{
    return textToFloat(progDef->getParam(param, subparam, arg));
}

float XmippProgram::getFloatParam(const char * param, int arg)
{
    return textToFloat(progDef->getParam(param, arg));
}

float XmippProgram::getFloatParam(const char * param, const char * subparam, int arg)
{
    return textToFloat(progDef->getParam(param, subparam, arg));
}

void XmippProgram::getListParam(const char * param, StringVector &list)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
    list.clear();
    for (size_t i = 0; i < paramDef->cmdArguments.size(); ++i)
        list.push_back(paramDef->cmdArguments[i]);
}

int XmippProgram::getCountParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
    return paramDef->cmdArguments.size();
}

bool XmippProgram::checkParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    if (paramDef == NULL)
        REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
    return paramDef->counter == 1;
}

bool XmippProgram::existsParam(const char * param)
{
    ParamDef * paramDef = progDef->findParam(param);
    return paramDef != NULL;
}


ParamDef * XmippProgram::getParamDef(const char * param) const
{
    return progDef->findParam(param);
}

const char * XmippProgram::name() const
{
    return progDef->name.c_str();
}

void XmippProgram::usage(int verb) const
{
    if (verbose)
    {
        ConsolePrinter cp;
        char * var = getenv("XMIPP_COLOR_OFF");
        if (var != NULL)
            cp.color = false;
        cp.printProgram(*progDef, verb);
    }
}

void XmippProgram::usage(const String & param, int verb)
{
    if (verbose)
    {
        ConsolePrinter cp;
        ParamDef * paramDef = progDef->findParam(param);
        if (paramDef == NULL)
            REPORT_ERROR(ERR_ARG_INCORRECT, ((String)"Doesn't exists param: " + param));
        cp.printParam(*paramDef, verb);
    }
    quit(0);
}

void XmippProgram::show() const
{}

int XmippProgram::version() const
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"");
}
