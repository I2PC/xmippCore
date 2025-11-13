/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#ifndef CORE_METADATA_STATIC_H
#define CORE_METADATA_STATIC_H

#include "metadata_label.h"
#include "metadata_row_sql.h"
#include "metadata_row_vec.h"

/** Static class to group some functions with labels.
 * This class holds several function to work with labels.
 * Also performs the static initialization of the labels data.
 */
class MDL
{
public:
    /** @name String conversions
     * @{
     */
    /** Converts an string to MDLabel */
    static void str2LabelVector(const String &labelsStr, std::vector<MDLabel> &labels);
    static MDLabel str2Label(const String &labelName);
    /** Converts MDLabel to string */
    static String label2Str(const MDLabel &label);
    /** Same as label2Str but escaping with '' to use in Sqlite. */
    static String label2StrSql(const MDLabel label);
    /** Converts MDLabel to string representing SQL column*/
    static String label2SqlColumn(const MDLabel label);
    /** Return the type of the label as String */
    static String labelType2Str(MDLabelType type);
    /** @} */

    /** @name Type checks
     * Several label type checks
     * @{
     */
    static bool isInt(const MDLabel label);
    static bool isLong(const MDLabel label);
    static bool isBool(const MDLabel label);
    static bool isString(const MDLabel label);
    static bool isDouble(const MDLabel label);
    static bool isVector(const MDLabel label);
    static bool isVectorLong(const MDLabel label);
    static bool isValidLabel(const MDLabel &label);
    static bool isValidLabel(const String &labelName);
    static MDLabelType labelType(const MDLabel label);
    static MDLabelType labelType(const String &labelName);
    static bool hasTag(const MDLabel label, const int tags);
    static bool isTextFile(const MDLabel label);
    static bool isMetadata(const MDLabel label);
    static bool isCtfParam(const MDLabel label);
    static bool isImage(const MDLabel label);
    static bool isVolume(const MDLabel label);
    static bool isStack(const MDLabel label);
    static bool isMicrograph(const MDLabel label);
    static bool isPSD(const MDLabel label);
    static std::map<String, MDLabel>& getLabelDict();
    /** @} */

    /** Add an alias for an existing label.
     * Params:
     *  label: The label id to create the alias to
     *  alias: The new alternative string name (alias)
     *  replace: If true, then new alias name will be used as the label
     *      string for writing back to file
     *  type: if provided, the label type will be replaced (only used with
     *      replace=True)
     *
     *  NOTE: Be aware that this function is not thread safe, so do not use
     *  it concurrently from different threads.
     */
    static void addLabelAlias(MDLabel label, const String &alias,
                              bool replace=false, MDLabelType type=LABEL_NOTYPE);

    /** Get an alias of an available BUFFER label with the provided name.
     * Params:
     *  alias: The name to create the alias from BUFFER label
     *  type: if provided, the label type will be replaced
     */
    static MDLabel getNewAlias(const String &alias, MDLabelType type=LABEL_NOTYPE);

    /** Reset the counter of buffer labels.
     * Use this function with care, since when you change the alias, the existing
     * metadata could change and write wrong values to file.
     */
    static void resetBufferIndex();

    static MDRowSql emptyHeaderSql();
    static MDRowVec emptyHeaderVec();
    static void emptifyHeader(MDRow&);

private:
    //Array of MDLabelData pointers
    static MDLabelData * data[MDL_LAST_LABEL+1];
    static std::map<std::string, MDLabel> names;
    static MDLabelStaticInit initialization; //Just for initialization
    static MDLabel bufferIndex; // Index that will be used to return next label alias

    /** Add predefined labels to be used in metadata */
    static void addLabel(MDLabel label, MDLabelType type, const String &name, int tags=TAGLABEL_NOTAG);
    /** This function will read extra label alias from XMIPP_LABEL_ALIASES var environment */
    static void addExtraAliases();
    friend class MDLabelStaticInit;

};//close class MLD definition


/** Just to work as static constructor for initialize labels data.
 */
class MDLabelStaticInit
{
private:
    MDLabelStaticInit()
    {
        //Just to be safe, initialize all data to null
        for (int i = (int)MDL_FIRST_LABEL; i <= (int)MDL_LAST_LABEL; ++i)
            MDL::data[i] = NULL;

        ///==== Add labels entries from here in the SAME ORDER as declared in ENUM ==========
        //The label MDL_OBJID is special and should not be used
        MDL::addLabel(MDL_OBJID, LABEL_SIZET, "objId");
        //The label MDL_GATHER_ID is special and should not be used
        MDL::addLabel(MDL_GATHER_ID, LABEL_SIZET, "gatherId");

        //MDL::addLabel(MDL_ANGLE_COMPARISON, LABEL_VECTOR_DOUBLE, "angle_comparison");
        //MDL::addLabelAlias(MDL_ANGLE_COMPARISON, "angleComparison"); //3.0

        MDL::addLabel(MDL_ANGLE_PSI, LABEL_DOUBLE, "anglePsi");
        MDL::addLabelAlias(MDL_ANGLE_PSI, "psi");
        MDL::addLabel(MDL_ANGLE_PSI2, LABEL_DOUBLE, "anglePsi2");
        MDL::addLabelAlias(MDL_ANGLE_PSI2, "psi2");
        MDL::addLabel(MDL_ANGLE_PSI3, LABEL_DOUBLE, "anglePsi3");
        MDL::addLabelAlias(MDL_ANGLE_PSI3, "psi3");
        MDL::addLabel(MDL_ANGLE_PSI_DIFF, LABEL_DOUBLE, "anglePsiDiff");
        MDL::addLabel(MDL_ANGLE_ROT, LABEL_DOUBLE, "angleRot");
        MDL::addLabelAlias(MDL_ANGLE_ROT, "rot");
        MDL::addLabel(MDL_ANGLE_ROT2, LABEL_DOUBLE, "angleRot2");
        MDL::addLabelAlias(MDL_ANGLE_ROT2, "rot2");
        MDL::addLabel(MDL_ANGLE_ROT3, LABEL_DOUBLE, "angleRot3");
        MDL::addLabelAlias(MDL_ANGLE_ROT3, "rot3");
        MDL::addLabel(MDL_ANGLE_ROT_DIFF, LABEL_DOUBLE, "angleRotDiff");
        MDL::addLabel(MDL_ANGLE_TILT, LABEL_DOUBLE, "angleTilt");
        MDL::addLabelAlias(MDL_ANGLE_TILT, "tilt");
        MDL::addLabel(MDL_ANGLE_TILT2, LABEL_DOUBLE, "angleTilt2");
        MDL::addLabelAlias(MDL_ANGLE_TILT2, "tilt2");
        MDL::addLabel(MDL_ANGLE_TILT3, LABEL_DOUBLE, "angleTilt3");
        MDL::addLabelAlias(MDL_ANGLE_TILT3, "tilt3");
        MDL::addLabel(MDL_ANGLE_TILT_DIFF, LABEL_DOUBLE, "angleTiltDiff");
        MDL::addLabel(MDL_ANGLE_DIFF0, LABEL_DOUBLE, "angleDiff0");
        MDL::addLabel(MDL_ANGLE_DIFF, LABEL_DOUBLE, "angleDiff");
        MDL::addLabel(MDL_ANGLE_DIFF2, LABEL_DOUBLE, "angleDiff2");
        MDL::addLabel(MDL_ANGLE_Y, LABEL_DOUBLE, "angleY");
        MDL::addLabel(MDL_ANGLE_Y2, LABEL_DOUBLE, "angleY2");
        MDL::addLabel(MDL_ANGLE_TEMPERATURE, LABEL_DOUBLE, "angleTemp");

        MDL::addLabel(MDL_APPLY_SHIFT, LABEL_BOOL, "applyShift");
        MDL::addLabel(MDL_AVG, LABEL_DOUBLE, "avg");
        MDL::addLabel(MDL_AVG_CHANGES_ORIENTATIONS, LABEL_DOUBLE, "avgChanOrient");
        MDL::addLabel(MDL_AVG_CHANGES_OFFSETS, LABEL_DOUBLE, "avgChanOffset");
        MDL::addLabel(MDL_AVG_CHANGES_CLASSES, LABEL_DOUBLE, "avgChanClass");
        MDL::addLabel(MDL_AVGPMAX, LABEL_DOUBLE, "avgPMax");

        MDL::addLabel(MDL_BFACTOR, LABEL_DOUBLE, "bFactor");
        MDL::addLabel(MDL_BGMEAN, LABEL_DOUBLE, "bgMean");
        MDL::addLabel(MDL_BLOCK_NUMBER, LABEL_INT, "blockNumber");

        MDL::addLabel(MDL_CL2D_CHANGES, LABEL_INT, "cl2dChanges");
        MDL::addLabel(MDL_CL2D_SIMILARITY, LABEL_DOUBLE, "cl2dSimilarity");
        MDL::addLabel(MDL_CLASS_COUNT, LABEL_SIZET, "classCount");
        MDL::addLabelAlias(MDL_CLASS_COUNT, "class_count"); //3.0
        MDL::addLabel(MDL_CLASS_PERCENTAGE, LABEL_DOUBLE, "classPercentage");
        MDL::addLabel(MDL_CLASSIFICATION_DATA, LABEL_VECTOR_DOUBLE, "classificationData");
        MDL::addLabelAlias(MDL_CLASSIFICATION_DATA, "ClassificationData");
        MDL::addLabel(MDL_CLASSIFICATION_DATA_SIZE, LABEL_SIZET, "classificationDatasize");
        MDL::addLabelAlias(MDL_CLASSIFICATION_DATA_SIZE, "ClassificationDataSize");
        MDL::addLabel(MDL_CLASSIFICATION_DPR_05, LABEL_DOUBLE, "classificationDPR05");
        MDL::addLabelAlias(MDL_CLASSIFICATION_DPR_05, "ClassificationDPR05");
        MDL::addLabel(MDL_CLASSIFICATION_FRC_05, LABEL_DOUBLE, "classificationFRC05");
        MDL::addLabelAlias(MDL_CLASSIFICATION_FRC_05, "ClassificationFRC05");
        MDL::addLabel(MDL_CLASSIFICATION_INTRACLASS_DISTANCE, LABEL_DOUBLE, "classificationIntraclassDistance");
        MDL::addLabelAlias(MDL_CLASSIFICATION_INTRACLASS_DISTANCE, "ClassificationIntraclassDistance");
        MDL::addLabel(MDL_COLOR, LABEL_INT, "color");
        MDL::addLabel(MDL_COMMENT, LABEL_STRING, "comment");
        MDL::addLabel(MDL_COST, LABEL_DOUBLE, "cost");
        MDL::addLabel(MDL_COST_PERCENTILE, LABEL_DOUBLE, "costPerc");
        MDL::addLabel(MDL_COORD_CONSENSUS_SCORE, LABEL_DOUBLE, "CoordConsScore");
        MDL::addLabel(MDL_COUNT2, LABEL_SIZET, "count2");
        MDL::addLabel(MDL_COUNT, LABEL_SIZET, "count");
        MDL::addLabel(MDL_CORR_DENOISED_PROJECTION, LABEL_DOUBLE, "corrDenoisedProjection");
        MDL::addLabel(MDL_CORR_DENOISED_NOISY, LABEL_DOUBLE, "corrDenoisedNoisy");
        MDL::addLabel(MDL_CORRELATION_IDX, LABEL_DOUBLE, "corrIdx");
        MDL::addLabel(MDL_CORRELATION_MASK, LABEL_DOUBLE, "corrMask");
        MDL::addLabel(MDL_CORRELATION_WEIGHT, LABEL_DOUBLE, "corrWeight");

        MDL::addLabel(MDL_CLASS_INTERSECTION_SIZE_PVALUE, LABEL_DOUBLE, "classIntersectionSizePValue");
        MDL::addLabel(MDL_CLASS_INTERSECTION_RELATIVE_SIZE_PVALUE, LABEL_DOUBLE, "classIntersectionRelativeSizePValue");

        MDL::addLabel(MDL_CRYSTAL_CELLX, LABEL_INT, "crystalCellx");
        MDL::addLabel(MDL_CRYSTAL_CELLY, LABEL_INT, "crystalCelly");
        MDL::addLabel(MDL_CRYSTAL_DISAPPEAR_THRE, LABEL_DOUBLE, "crystalDisthresh");
        MDL::addLabel(MDL_CRYSTAL_LATTICE_A, LABEL_VECTOR_DOUBLE, "crystalLatticeA");
        MDL::addLabel(MDL_CRYSTAL_LATTICE_B, LABEL_VECTOR_DOUBLE, "crystalLatticeB");
        MDL::addLabel(MDL_CRYSTAL_ORTHO_PRJ, LABEL_BOOL, "crystalOrthoProj");
        MDL::addLabel(MDL_CRYSTAL_PROJ, LABEL_BOOL, "crystalProj");
        MDL::addLabel(MDL_CRYSTAL_SHFILE, LABEL_STRING, "crystalShiftFile");
        MDL::addLabel(MDL_CRYSTAL_SHIFTX, LABEL_DOUBLE, "crystalShiftX");
        MDL::addLabel(MDL_CRYSTAL_SHIFTY, LABEL_DOUBLE, "crystalShiftY");
        MDL::addLabel(MDL_CRYSTAL_SHIFTZ, LABEL_DOUBLE, "crystalShiftZ");
        MDL::addLabel(MDL_CRYSTAL_NOISE_SHIFT,LABEL_VECTOR_DOUBLE, "crystalNoiseShift");
        MDL::addLabel(MDL_CTF_BG_BASELINE, LABEL_DOUBLE, "ctfBgBaseline");
        MDL::addLabelAlias(MDL_CTF_BG_BASELINE, "CTFBG_Baseline");//3.0
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_ANGLE, LABEL_DOUBLE, "ctfBgGaussian2Angle");
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_ANGLE, "CTFBG_Gaussian2_Angle"); //3.0
        MDL::addLabel(MDL_CTF_BG_R1, LABEL_DOUBLE, "ctfBgR1");
        MDL::addLabel(MDL_CTF_BG_R2, LABEL_DOUBLE, "ctfBgR2");
        MDL::addLabel(MDL_CTF_BG_R3, LABEL_DOUBLE, "ctfBgR3");

        MDL::addLabel(MDL_CTF_DATA_PHASE_FLIPPED, LABEL_BOOL, "ctfPhaseFlipped");
        MDL::addLabel(MDL_CTF_CORRECTED, LABEL_BOOL, "ctfCorrected");
        MDL::addLabel(MDL_CTF_X0, LABEL_DOUBLE, "ctfX0");
        MDL::addLabel(MDL_CTF_XF, LABEL_DOUBLE, "ctfXF");
        MDL::addLabel(MDL_CTF_Y0, LABEL_DOUBLE, "ctfY0");
        MDL::addLabel(MDL_CTF_YF, LABEL_DOUBLE, "ctfYF");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEUA, LABEL_DOUBLE, "ctfDefocusPlaneUA");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEUB, LABEL_DOUBLE, "ctfDefocusPlaneUB");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEUC, LABEL_DOUBLE, "ctfDefocusPlaneUC");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEVA, LABEL_DOUBLE, "ctfDefocusPlaneVA");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEVB, LABEL_DOUBLE, "ctfDefocusPlaneVB");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEVC, LABEL_DOUBLE, "ctfDefocusPlaneVC");


        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_CU, LABEL_DOUBLE, "ctfBgGaussian2CU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_CV, LABEL_DOUBLE, "ctfBgGaussian2CV");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_K, LABEL_DOUBLE, "ctfBgGaussian2K");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAU, LABEL_DOUBLE, "ctfBgGaussian2SigmaU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAV, LABEL_DOUBLE, "ctfBgGaussian2SigmaV");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_ANGLE, LABEL_DOUBLE, "ctfBgGaussianAngle");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_CU, LABEL_DOUBLE, "ctfBgGaussianCU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_CV, LABEL_DOUBLE, "ctfBgGaussianCV");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_K, LABEL_DOUBLE, "ctfBgGaussianK");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_SIGMAU, LABEL_DOUBLE, "ctfBgGaussianSigmaU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_SIGMAV, LABEL_DOUBLE, "ctfBgGaussianSigmaV");
        MDL::addLabel(MDL_CTF_PHASE_SHIFT, LABEL_DOUBLE, "ctfVPPphaseshift");
        MDL::addLabel(MDL_CTF_VPP_RADIUS, LABEL_DOUBLE, "ctfVPPRadius");

        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_CU, "CTFBG_Gaussian2_CU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_CV, "CTFBG_Gaussian2_CV");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_K, "CTFBG_Gaussian2_K");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_SIGMAU, "CTFBG_Gaussian2_SigmaU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_SIGMAV, "CTFBG_Gaussian2_SigmaV");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_ANGLE, "CTFBG_Gaussian_Angle");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_CU, "CTFBG_Gaussian_CU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_CV, "CTFBG_Gaussian_CV");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_K, "CTFBG_Gaussian_K");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_SIGMAU, "CTFBG_Gaussian_SigmaU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_SIGMAV, "CTFBG_Gaussian_SigmaV");//3.0

        MDL::addLabel(MDL_CTF_BG_SQRT_ANGLE, LABEL_DOUBLE, "ctfBgSqrtAngle");
        MDL::addLabel(MDL_CTF_BG_SQRT_K, LABEL_DOUBLE, "ctfBgSqrtK");
        MDL::addLabel(MDL_CTF_BG_SQRT_U, LABEL_DOUBLE, "ctfBgSqrtU");
        MDL::addLabel(MDL_CTF_BG_SQRT_V, LABEL_DOUBLE, "ctfBgSqrtV");
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_ANGLE, "CTFBG_Sqrt_Angle");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_K, "CTFBG_Sqrt_K");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_U, "CTFBG_Sqrt_U");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_V, "CTFBG_Sqrt_V");  //3.0

        MDL::addLabel(MDL_CONTINUOUS_X, LABEL_DOUBLE, "continuousX");
        MDL::addLabel(MDL_CONTINUOUS_Y, LABEL_DOUBLE, "continuousY");
        MDL::addLabel(MDL_CONTINUOUS_FLIP, LABEL_BOOL, "continuousFlip");
        MDL::addLabel(MDL_CONTINUOUS_GRAY_A, LABEL_DOUBLE, "continuousA");
        MDL::addLabel(MDL_CONTINUOUS_GRAY_B, LABEL_DOUBLE, "continuousB");
        MDL::addLabel(MDL_CONTINUOUS_SCALE_ANGLE, LABEL_DOUBLE, "continuousScaleAngle");
        MDL::addLabel(MDL_CONTINUOUS_SCALE_X, LABEL_DOUBLE, "continuousScaleX");
        MDL::addLabel(MDL_CONTINUOUS_SCALE_Y, LABEL_DOUBLE, "continuousScaleY");
        MDL::addLabel(MDL_CTF_CA, LABEL_DOUBLE, "ctfChromaticAberration");
        MDL::addLabel(MDL_CTF_CONVERGENCE_CONE, LABEL_DOUBLE, "ctfConvergenceCone");
        MDL::addLabel(MDL_CTF_CRIT_NONASTIGMATICVALIDITY, LABEL_DOUBLE, "ctfCritNonAstigmaticValidty");
        MDL::addLabel(MDL_CTF_CRIT_DAMPING, LABEL_DOUBLE, "ctfCritDamping");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZEROAVG, LABEL_DOUBLE, "ctfCritFirstZero");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, LABEL_DOUBLE, "ctfCritDisagree");
        MDL::addLabel(MDL_CTF_CRIT_MAXFREQ, LABEL_DOUBLE, "ctfCritMaxFreq");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZERORATIO, LABEL_DOUBLE, "ctfCritfirstZeroRatio");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_RATIO, LABEL_DOUBLE, "ctfCritFirstMinFirstZeroRatio");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO, LABEL_DOUBLE, "ctfCritCtfMargin");
        MDL::addLabel(MDL_CTF_CRIT_FITTINGCORR13, LABEL_DOUBLE, "ctfCritCorr13");
        MDL::addLabel(MDL_CTF_CRIT_FITTINGSCORE, LABEL_DOUBLE, "ctfCritFitting");
        MDL::addLabel(MDL_CTF_CRIT_ICENESS, LABEL_DOUBLE, "ctfCritIceness");
        MDL::addLabel(MDL_CTF_CRIT_NORMALITY, LABEL_DOUBLE, "ctfCritNormality");
        MDL::addLabel(MDL_CTF_CRIT_PSDCORRELATION90, LABEL_DOUBLE, "ctfCritPsdCorr90");
        MDL::addLabel(MDL_CTF_CRIT_PSDPCA1VARIANCE, LABEL_DOUBLE, "ctfCritPsdPCA1");
        MDL::addLabel(MDL_CTF_CRIT_PSDPCARUNSTEST, LABEL_DOUBLE, "ctfCritPsdPCARuns");
        MDL::addLabel(MDL_CTF_CRIT_PSDRADIALINTEGRAL, LABEL_DOUBLE, "ctfCritPsdInt");
        MDL::addLabel(MDL_CTF_CRIT_PSDVARIANCE, LABEL_DOUBLE, "ctfCritPsdStdQ");
        MDL::addLabel(MDL_CTF_CS, LABEL_DOUBLE, "ctfSphericalAberration");
        MDL::addLabel(MDL_CTF_DEFOCUSA, LABEL_DOUBLE, "ctfDefocusA");//average defocus
        MDL::addLabel(MDL_CTF_DEFOCUS_ANGLE, LABEL_DOUBLE, "ctfDefocusAngle");
        MDL::addLabel(MDL_CTF_DEFOCUSU, LABEL_DOUBLE, "ctfDefocusU");
        MDL::addLabel(MDL_CTF_DEFOCUSV, LABEL_DOUBLE, "ctfDefocusV");
        MDL::addLabel(MDL_CTF_DEFOCUS_CHANGE, LABEL_DOUBLE, "ctfDefocusChange");
        MDL::addLabel(MDL_CTF_DEFOCUS_R2, LABEL_DOUBLE, "ctfDefocusR2");
        MDL::addLabel(MDL_CTF_DEFOCUS_RESIDUAL, LABEL_DOUBLE, "ctfDefocusResidual");
        MDL::addLabel(MDL_CTF_DEFOCUS_COEFS, LABEL_VECTOR_DOUBLE, "ctfDefocusCoeficients");
        MDL::addLabel(MDL_CTF_DIMENSIONS, LABEL_VECTOR_DOUBLE, "ctfDimensions");
        MDL::addLabel(MDL_CTF_DOWNSAMPLE_PERFORMED, LABEL_DOUBLE, "CtfDownsampleFactor");
        MDL::addLabel(MDL_CTF_ENERGY_LOSS, LABEL_DOUBLE, "ctfEnergyLoss");
        MDL::addLabel(MDL_CTF_ENVELOPE, LABEL_DOUBLE, "ctfEnvelope");
        MDL::addLabel(MDL_CTF_ENVELOPE_PLOT, LABEL_STRING, "ctfEnvelopePlot");
        MDL::addLabel(MDL_CTF_GROUP, LABEL_INT, "ctfGroup");
        MDL::addLabel(MDL_CTF_INPUTPARAMS, LABEL_STRING, "ctfInputParams", TAGLABEL_TEXTFILE);
        MDL::addLabel(MDL_CTF_K, LABEL_DOUBLE, "ctfK");
        MDL::addLabel(MDL_CTF_ENV_R0, LABEL_DOUBLE, "ctfEnvR0");
        MDL::addLabel(MDL_CTF_ENV_R1, LABEL_DOUBLE, "ctfEnvR1");
        MDL::addLabel(MDL_CTF_ENV_R2, LABEL_DOUBLE, "ctfEnvR2");
        MDL::addLabel(MDL_CTF_LAMBDA, LABEL_DOUBLE, "ctfLambda");
        MDL::addLabel(MDL_CTF_LENS_STABILITY, LABEL_DOUBLE, "ctfLensStability");
        MDL::addLabel(MDL_CTF_LONGITUDINAL_DISPLACEMENT, LABEL_DOUBLE, "ctfLongitudinalDisplacement");
        MDL::addLabel(MDL_CTF_MODEL2, LABEL_STRING, "ctfModel2", TAGLABEL_CTFPARAM);
        MDL::addLabel(MDL_CTF_MODEL, LABEL_STRING, "ctfModel", TAGLABEL_CTFPARAM);
        MDL::addLabel(MDL_CTF_Q0, LABEL_DOUBLE, "ctfQ0");
        MDL::addLabel(MDL_CTF_SAMPLING_RATE, LABEL_DOUBLE, "ctfSamplingRate");
        MDL::addLabel(MDL_CTF_SAMPLING_RATE_Z, LABEL_DOUBLE, "ctfSamplingRateZ");
        MDL::addLabel(MDL_CTF_TRANSVERSAL_DISPLACEMENT, LABEL_DOUBLE, "ctfTransversalDisplacement");
        MDL::addLabel(MDL_CTF_VOLTAGE, LABEL_DOUBLE, "ctfVoltage");
        MDL::addLabel(MDL_CTF_XRAY_LENS_TYPE, LABEL_STRING, "ctfXrayLensType");
        MDL::addLabel(MDL_CTF_XRAY_OUTER_ZONE_WIDTH, LABEL_DOUBLE, "ctfXrayOuterZoneWidth");
        MDL::addLabel(MDL_CTF_XRAY_ZONES_NUMBER, LABEL_DOUBLE, "ctfXrayZonesN");
        MDL::addLabel(MDL_CUMULATIVE_SSNR, LABEL_DOUBLE, "cumulativeSSNR");

        MDL::addLabelAlias(MDL_CTF_CA, "CTF_Chromatic_aberration"); //3.0
        MDL::addLabelAlias(MDL_CTF_CONVERGENCE_CONE, "CTF_Convergence_cone"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_DAMPING, "CTFCrit_damping"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FIRSTZEROAVG, "CTFCrit_firstZero"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, "CTFCrit_disagree"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FIRSTZERORATIO, "CTFCrit_firstZeroRatio"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FITTINGCORR13, "CTFCrit_Corr13"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FITTINGSCORE, "CTFCrit_Fitting"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_NORMALITY, "CTFCrit_Normality"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDCORRELATION90, "CTFCrit_psdcorr90"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDPCA1VARIANCE, "CTFCrit_PSDPCA1"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDPCARUNSTEST, "CTFCrit_PSDPCARuns"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDRADIALINTEGRAL, "CTFCrit_psdint"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDVARIANCE, "CTFCrit_PSDStdQ"); //3.0
        MDL::addLabelAlias(MDL_CTF_CS, "CTF_Spherical_aberration"); //3.0
        MDL::addLabelAlias(MDL_CTF_DEFOCUSA, "CTF_Defocus_A"); //3.0//average defocus
        MDL::addLabelAlias(MDL_CTF_DEFOCUS_ANGLE, "CTF_Defocus_angle"); //3.0
        MDL::addLabelAlias(MDL_CTF_DEFOCUSU, "CTF_Defocus_U"); //3.0
        MDL::addLabelAlias(MDL_CTF_DEFOCUSV, "CTF_Defocus_V"); //3.0
        MDL::addLabelAlias(MDL_CTF_DIMENSIONS, "CTF_Xray_dimensions"); //3.0
        MDL::addLabelAlias(MDL_CTF_DOWNSAMPLE_PERFORMED, "CTFDownsampleFactor"); //3.0
        MDL::addLabelAlias(MDL_CTF_ENERGY_LOSS, "CTF_Energy_loss"); //3.0
        MDL::addLabelAlias(MDL_CTF_GROUP, "CTFGroup"); //3.0
        MDL::addLabelAlias(MDL_CTF_INPUTPARAMS, "CTFInputParams"); //3.0
        MDL::addLabelAlias(MDL_CTF_K, "CTF_K"); //3.0
        MDL::addLabelAlias(MDL_CTF_LAMBDA, "CTF_Xray_lambda"); //3.0
        MDL::addLabelAlias(MDL_CTF_LENS_STABILITY, "CTF_Lens_stability"); //3.0
        MDL::addLabelAlias(MDL_CTF_LONGITUDINAL_DISPLACEMENT, "CTF_Longitudinal_displacement"); //3.0
        MDL::addLabelAlias(MDL_CTF_MODEL2, "CTFModel2"); //3.0
        MDL::addLabelAlias(MDL_CTF_MODEL, "CTFModel"); //3.0
        MDL::addLabelAlias(MDL_CTF_Q0, "CTF_Q0"); //3.0
        MDL::addLabelAlias(MDL_CTF_SAMPLING_RATE, "CTF_Sampling_rate"); //3.0
        MDL::addLabelAlias(MDL_CTF_SAMPLING_RATE_Z, "CTF_Sampling_rate_z"); //3.0
        MDL::addLabelAlias(MDL_CTF_TRANSVERSAL_DISPLACEMENT, "CTF_Transversal_displacement"); //3.0
        MDL::addLabelAlias(MDL_CTF_VOLTAGE, "CTF_Voltage"); //3.0
        MDL::addLabelAlias(MDL_CTF_XRAY_LENS_TYPE, "CTF_Xray_lens_type"); //3.0
        MDL::addLabelAlias(MDL_CTF_XRAY_OUTER_ZONE_WIDTH, "CTF_Xray_OuterZoneWidth"); //3.0
        MDL::addLabelAlias(MDL_CTF_XRAY_ZONES_NUMBER, "CTF_Xray_ZonesN"); //3.0

        MDL::addLabel(MDL_DATATYPE, LABEL_INT, "datatype");
        MDL::addLabel(MDL_DEFGROUP, LABEL_INT, "defocusGroup");
        MDL::addLabel(MDL_DIMENSIONS_2D, LABEL_VECTOR_DOUBLE, "dimensions2D");
        MDL::addLabel(MDL_DIMENSIONS_3D, LABEL_VECTOR_DOUBLE, "dimensions3D");
        MDL::addLabel(MDL_DIMRED, LABEL_VECTOR_DOUBLE, "dimredCoeffs");
        MDL::addLabel(MDL_DIRECTION, LABEL_VECTOR_DOUBLE, "direction");
        MDL::addLabel(MDL_DM3_IDTAG, LABEL_INT, "dm3IdTag");
        MDL::addLabel(MDL_DM3_NODEID, LABEL_INT, "dm3NodeId");
        MDL::addLabel(MDL_DM3_NUMBER_TYPE, LABEL_INT, "dm3NumberType");
        MDL::addLabel(MDL_DM3_PARENTID, LABEL_INT, "dm3ParentID");
        MDL::addLabel(MDL_DM3_SIZE, LABEL_INT, "dm3Size");
        MDL::addLabel(MDL_DM3_TAGCLASS, LABEL_STRING, "dm3TagClass");
        MDL::addLabel(MDL_DM3_TAGNAME, LABEL_STRING, "dm3TagName");
        MDL::addLabel(MDL_DM3_VALUE, LABEL_VECTOR_DOUBLE, "dm3Value");
        MDL::addLabel(MDL_DOSE, LABEL_DOUBLE, "dose");

        MDL::addLabel(MDL_ENABLED, LABEL_INT, "enabled");

        //MDL_EXECUTION_DATE so far an string but may change...
        MDL::addLabel(MDL_DATE, LABEL_STRING, "date");
        MDL::addLabel(MDL_TIME, LABEL_DOUBLE, "time");

        MDL::addLabel(MDL_FLIP, LABEL_BOOL, "flip");
        MDL::addLabelAlias(MDL_FLIP, "Flip");
        MDL::addLabel(MDL_FOM, LABEL_DOUBLE, "fom");
        MDL::addLabel(MDL_FRAME_ID, LABEL_SIZET, "frameId");

        MDL::addLabel(MDL_GRAPH_DISTANCE2MAX, LABEL_DOUBLE, "distance2MaxGraph");
        MDL::addLabel(MDL_GRAPH_DISTANCE2MAX_PREVIOUS, LABEL_DOUBLE, "distance2MaxGraphPrevious");
        MDL::addLabel(MDL_GRAPH_CC, LABEL_DOUBLE, "graphCC");
        MDL::addLabel(MDL_GRAPH_CC_PREVIOUS, LABEL_DOUBLE, "graphCCPrevious");
        MDL::addLabel(MDL_HALF1, LABEL_STRING, "half1");
        MDL::addLabel(MDL_HALF2, LABEL_STRING, "half2");
        MDL::addLabel(MDL_ASSIGNED_DIR_REF_CC, LABEL_DOUBLE, "assignedDirRefCC");

        MDL::addLabel(MDL_IDX, LABEL_SIZET, "index");
        MDL::addLabel(MDL_IMAGE1, LABEL_STRING, "image1", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE2, LABEL_STRING, "image2", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE3, LABEL_STRING, "image3", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE4, LABEL_STRING, "image4", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE5, LABEL_STRING, "image5", TAGLABEL_IMAGE);

        MDL::addLabelAlias(MDL_IMAGE1, "associatedImage1"); //3.0
        MDL::addLabelAlias(MDL_IMAGE2, "associatedImage2"); //3.0
        MDL::addLabelAlias(MDL_IMAGE3, "associatedImage3"); //3.0
        MDL::addLabelAlias(MDL_IMAGE4, "associatedImage4"); //3.0
        MDL::addLabelAlias(MDL_IMAGE5, "associatedImage5"); //3.0

        MDL::addLabel(MDL_IMAGE, LABEL_STRING, "image", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_COVARIANCE, LABEL_STRING, "imageCovariance", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_IDX, LABEL_SIZET, "imageIndex");
        MDL::addLabel(MDL_IMAGE_ORIGINAL, LABEL_STRING, "imageOriginal", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_REF, LABEL_STRING, "imageRef", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_RESIDUAL, LABEL_STRING, "imageResidual", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_TILTED, LABEL_STRING, "imageTilted", TAGLABEL_IMAGE);

        MDL::addLabelAlias(MDL_IMAGE_ORIGINAL, "original_image"); //3.0
        MDL::addLabelAlias(MDL_IMAGE_TILTED, "tilted_image"); //3.0

        MDL::addLabel(MDL_IMED, LABEL_DOUBLE, "imedValue");

        MDL::addLabel(MDL_IMGMD, LABEL_STRING, "imageMetaData", TAGLABEL_METADATA);
        MDL::addLabel(MDL_INTSCALE, LABEL_DOUBLE, "intScale");

        MDL::addLabel(MDL_ITEM_ID, LABEL_SIZET, "itemId");
        MDL::addLabel(MDL_ITER, LABEL_INT, "iterationNumber");

        MDL::addLabel(MDL_KERDENSOM_FUNCTIONAL, LABEL_DOUBLE, "kerdensomFunctional");
        MDL::addLabel(MDL_KERDENSOM_REGULARIZATION, LABEL_DOUBLE, "kerdensomRegularization");
        MDL::addLabel(MDL_KERDENSOM_SIGMA, LABEL_DOUBLE, "kerdensomSigma");

        MDL::addLabel(MDL_KEYWORDS, LABEL_STRING, "keywords");
        MDL::addLabel(MDL_KSTEST, LABEL_DOUBLE, "kstest");
        MDL::addLabel(MDL_LL, LABEL_DOUBLE, "logLikelihood");
        MDL::addLabelAlias(MDL_LL, "LL");
        MDL::addLabel(MDL_LOCAL_ALIGNMENT_PATCHES, LABEL_VECTOR_SIZET, "localAlignmentPatches");
        MDL::addLabel(MDL_LOCAL_ALIGNMENT_COEFFS_X, LABEL_VECTOR_DOUBLE, "localAlignmentCoeffsX");
        MDL::addLabel(MDL_LOCAL_ALIGNMENT_COEFFS_Y, LABEL_VECTOR_DOUBLE, "localAlignmentCoeffsY");
        MDL::addLabel(MDL_LOCAL_ALIGNMENT_CONF_2_5_PERC, LABEL_DOUBLE, "localAlignmnentConf25Perc");
        MDL::addLabel(MDL_LOCAL_ALIGNMENT_CONF_97_5_PERC, LABEL_DOUBLE, "localAlignmnentConf955Perc");
        MDL::addLabel(MDL_LOCAL_ALIGNMENT_CONTROL_POINTS, LABEL_VECTOR_SIZET, "localAlignmentControlPoints");
        MDL::addLabel(MDL_LOCAL_AVERAGE, LABEL_DOUBLE, "localAverage");
        MDL::addLabel(MDL_MACRO_CMD, LABEL_STRING, "macroCmd");
        MDL::addLabel(MDL_MACRO_CMD_ARGS, LABEL_STRING, "macroCmdArgs");
        MDL::addLabel(MDL_MAGNIFICATION, LABEL_DOUBLE, "magnification");
        MDL::addLabel(MDL_MAPTOPOLOGY, LABEL_STRING, "mapTopology");
        MDL::addLabel(MDL_MASK, LABEL_STRING, "mask", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_MAXCC, LABEL_DOUBLE, "maxCC");
        MDL::addLabel(MDL_MAXCC_PERCENTILE, LABEL_DOUBLE, "maxCCPerc");
        MDL::addLabel(MDL_MAX, LABEL_DOUBLE, "max");
        MDL::addLabel(MDL_MAXCC_PREVIOUS, LABEL_DOUBLE, "maxCCprevious");
        MDL::addLabel(MDL_MICROGRAPH_ID, LABEL_SIZET, "micrographId");
        MDL::addLabel(MDL_MICROGRAPH, LABEL_STRING, "micrograph", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_MOVIE_ID, LABEL_SIZET, "micrographMovieId");
        MDL::addLabel(MDL_MICROGRAPH_MOVIE, LABEL_STRING, "movie", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_PARTICLES, LABEL_STRING, "micrographParticles", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_ORIGINAL, LABEL_STRING, "micrographOriginal", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_TILTED, LABEL_STRING, "micrographTilted", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_TILTED_ORIGINAL, LABEL_STRING, "micrographTiltedOriginal", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MIN, LABEL_DOUBLE, "min");
        MDL::addLabel(MDL_MIRRORFRAC, LABEL_DOUBLE, "mirrorFraction");
        MDL::addLabel(MDL_MISSINGREGION_NR, LABEL_INT, "missingRegionNumber");
        MDL::addLabelAlias(MDL_MISSINGREGION_NR, "Wedge");
        MDL::addLabel(MDL_MISSINGREGION_THX0, LABEL_DOUBLE, "missingRegionThetaX0");
        MDL::addLabel(MDL_MISSINGREGION_THXF, LABEL_DOUBLE, "missingRegionThetaXF");

        MDL::addLabel(MDL_MLF_CTF,    LABEL_DOUBLE, "mlfCtf");
        MDL::addLabel(MDL_MLF_WIENER, LABEL_DOUBLE, "mlfWiener");
        MDL::addLabel(MDL_MLF_SIGNAL, LABEL_DOUBLE, "mlfSignal");
        MDL::addLabel(MDL_MLF_NOISE,  LABEL_DOUBLE, "mlfNoise");

        MDL::addLabel(MDL_MISSINGREGION_THY0, LABEL_DOUBLE, "missingRegionThetaY0");
        MDL::addLabel(MDL_MISSINGREGION_THYF, LABEL_DOUBLE, "missingRegionThetaYF");
        MDL::addLabel(MDL_MISSINGREGION_TYPE, LABEL_STRING, "missingRegionType");
        MDL::addLabel(MDL_MODELFRAC, LABEL_DOUBLE, "modelFraction");
        MDL::addLabel(MDL_NEIGHBORHOOD_RADIUS, LABEL_DOUBLE, "neighborhoodRadius");
        MDL::addLabel(MDL_NEIGHBOR, LABEL_SIZET, "neighbor");
        MDL::addLabel(MDL_NEIGHBORS, LABEL_VECTOR_SIZET, "neighbors");
        MDL::addLabel(MDL_NMA, LABEL_VECTOR_DOUBLE, "nmaDisplacements");
        MDL::addLabelAlias(MDL_NMA, "NMADisplacements");//3.0
        MDL::addLabel(MDL_NMA_ATOMSHIFT, LABEL_DOUBLE, "nmaAtomShift");
        MDL::addLabel(MDL_NMA_COLLECTIVITY, LABEL_DOUBLE, "nmaCollectivity");
        MDL::addLabel(MDL_NMA_ENERGY, LABEL_DOUBLE, "nmaEnergy");
        MDL::addLabel(MDL_NMA_MINRANGE, LABEL_DOUBLE, "nmaMin");
        MDL::addLabel(MDL_NMA_MAXRANGE, LABEL_DOUBLE, "nmaMax");
        MDL::addLabel(MDL_NMA_MODEFILE, LABEL_STRING, "nmaModefile", TAGLABEL_TEXTFILE);
        MDL::addLabelAlias(MDL_NMA_MODEFILE, "NMAModefile");//3.0
        MDL::addLabel(MDL_NMA_SCORE, LABEL_DOUBLE, "nmaScore");
        MDL::addLabel(MDL_NMA_EIGENVAL, LABEL_DOUBLE, "nmaEigenval");
        MDL::addLabel(MDL_NOISE_ANGLES, LABEL_VECTOR_DOUBLE, "noiseAngles");
        MDL::addLabel(MDL_NOISE_COORD, LABEL_VECTOR_DOUBLE, "noiseCoord");
        MDL::addLabel(MDL_NOISE_PARTICLE_COORD, LABEL_VECTOR_DOUBLE, "noiseParticleCoord");
        MDL::addLabel(MDL_NOISE_PIXEL_LEVEL, LABEL_VECTOR_DOUBLE, "noisePixelLevel");
        MDL::addLabel(MDL_ORDER, LABEL_SIZET, "order_");
        MDL::addLabel(MDL_ORIGIN_X, LABEL_DOUBLE, "originX");
        MDL::addLabel(MDL_ORIGIN_Y, LABEL_DOUBLE, "originY");
        MDL::addLabel(MDL_ORIGIN_Z, LABEL_DOUBLE, "originZ");
        MDL::addLabel(MDL_PARTICLE_ID, LABEL_SIZET, "particleId");
        MDL::addLabel(MDL_PHANTOM_BGDENSITY, LABEL_DOUBLE, "phantomBGDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_CENTER, LABEL_VECTOR_DOUBLE, "featureCenter");
        MDL::addLabel(MDL_PHANTOM_FEATURE_DENSITY, LABEL_DOUBLE, "featureDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_OPERATION, LABEL_STRING, "featureOperation");
        MDL::addLabel(MDL_PHANTOM_FEATURE_SPECIFIC, LABEL_VECTOR_DOUBLE, "featureSpecificVector");
        MDL::addLabel(MDL_PHANTOM_FEATURE_TYPE, LABEL_STRING, "featureType");
        MDL::addLabel(MDL_PHANTOM_SCALE, LABEL_DOUBLE, "phantomScale");

        MDL::addLabel(MDL_OPTICALFLOW_MEANX, LABEL_DOUBLE, "opticalMeanX");
        MDL::addLabel(MDL_OPTICALFLOW_MEANY, LABEL_DOUBLE, "opticalMeanY");
        MDL::addLabel(MDL_OPTICALFLOW_STDX, LABEL_DOUBLE, "opticalStdX");
        MDL::addLabel(MDL_OPTICALFLOW_STDY, LABEL_DOUBLE, "opticalStdY");

        MDL::addLabel(MDL_PICKING_STATE, LABEL_STRING, "pickingState");
        MDL::addLabelAlias(MDL_PICKING_STATE, "picking_state");//3.0
        MDL::addLabel(MDL_PICKING_MICROGRAPH_STATE, LABEL_STRING, "pickingMicrographState");
        MDL::addLabelAlias(MDL_PICKING_MICROGRAPH_STATE, "micrograph_state");//3.0
        MDL::addLabel(MDL_PICKING_PARTICLE_SIZE, LABEL_INT, "particleSize");
        MDL::addLabel(MDL_PICKING_AUTOPICKPERCENT, LABEL_INT, "autopickPercent");
        MDL::addLabel(MDL_PICKING_TEMPLATES, LABEL_INT, "templatesNum");
        MDL::addLabel(MDL_PICKING_AUTOPARTICLES_SIZE, LABEL_INT, "autoParticlesNum");
        MDL::addLabel(MDL_PICKING_MANUALPARTICLES_SIZE, LABEL_INT, "manualParticlesNum");

        MDL::addLabel(MDL_PMAX, LABEL_DOUBLE, "pMax");
        MDL::addLabelAlias(MDL_PMAX, "Pmax");
        MDL::addLabelAlias(MDL_PMAX, "sumP");
        MDL::addLabel(MDL_POINTSASYMETRICUNIT, LABEL_SIZET, "pointsAsymmetricUnit");
        MDL::addLabelAlias(MDL_POINTSASYMETRICUNIT, "pointsasymmetricUnit");
        MDL::addLabel(MDL_PRJ_ANGFILE, LABEL_STRING, "projAngleFile");
        MDL::addLabelAlias(MDL_PRJ_ANGFILE, "angleFile");//3.0
        MDL::addLabel(MDL_PRJ_DIMENSIONS, LABEL_VECTOR_DOUBLE, "projDimensions");
        MDL::addLabel(MDL_PRJ_PSI_NOISE, LABEL_VECTOR_DOUBLE, "projPsiNoise");
        MDL::addLabel(MDL_PRJ_PSI_RANDSTR, LABEL_STRING, "projPsiRandomness");
        MDL::addLabel(MDL_PRJ_PSI_RANGE, LABEL_VECTOR_DOUBLE, "projPsiRange");
        MDL::addLabel(MDL_PRJ_ROT_NOISE, LABEL_VECTOR_DOUBLE, "projRotNoise");
        MDL::addLabel(MDL_PRJ_ROT_RANDSTR, LABEL_STRING, "projRotRandomness");
        MDL::addLabel(MDL_PRJ_ROT_RANGE, LABEL_VECTOR_DOUBLE, "projRotRange");
        MDL::addLabel(MDL_PRJ_TILT_NOISE, LABEL_VECTOR_DOUBLE, "projTiltNoise");
        MDL::addLabel(MDL_PRJ_TILT_RANDSTR, LABEL_STRING, "projTiltRandomness");
        MDL::addLabel(MDL_PRJ_TILT_RANGE, LABEL_VECTOR_DOUBLE, "projTiltRange");
        MDL::addLabel(MDL_PRJ_VOL, LABEL_STRING, "projVolume", TAGLABEL_VOLUME);

        MDL::addLabel(MDL_PROGRAM, LABEL_STRING, "program");
        MDL::addLabel(MDL_USER, LABEL_STRING, "user");

        MDL::addLabel(MDL_PSD_ENHANCED, LABEL_STRING, "psdEnhanced", TAGLABEL_IMAGE);
        MDL::addLabelAlias(MDL_PSD_ENHANCED, "enhancedPowerSpectrum");//3.0
        MDL::addLabel(MDL_PSD, LABEL_STRING, "psd", TAGLABEL_PSD);
        MDL::addLabelAlias(MDL_PSD, "powerSpectrum");//3.0
        MDL::addLabel(MDL_RANDOMSEED, LABEL_INT, "randomSeed");
        MDL::addLabel(MDL_REF2, LABEL_INT, "ref2");
        MDL::addLabel(MDL_REF3D, LABEL_INT, "ref3d");
        MDL::addLabel(MDL_REF, LABEL_INT, "ref");
        MDL::addLabelAlias(MDL_REF, "Ref");
        MDL::addLabel(MDL_REFMD, LABEL_STRING, "referenceMetaData", TAGLABEL_METADATA);

        MDL::addLabel(MDL_RESIDUE, LABEL_INT, "residue");
        MDL::addLabel(MDL_RESOLUTION_ANISOTROPY, LABEL_DOUBLE, "resolutionAnisotropy");
        MDL::addLabel(MDL_RESOLUTION_DPR, LABEL_DOUBLE, "resolutionDPR");
        MDL::addLabel(MDL_RESOLUTION_ERRORL2, LABEL_DOUBLE, "resolutionErrorL2");
        MDL::addLabel(MDL_RESOLUTION_FSO, LABEL_DOUBLE, "resolutionFSO");
        MDL::addLabel(MDL_RESOLUTION_FRC, LABEL_DOUBLE, "resolutionFRC");
        MDL::addLabel(MDL_RESOLUTION_FRCRANDOMNOISE, LABEL_DOUBLE, "resolutionFRCRandomNoise");
        MDL::addLabel(MDL_RESOLUTION_FREQ, LABEL_DOUBLE, "resolutionFreqFourier");
        MDL::addLabel(MDL_RESOLUTION_FREQ2, LABEL_DOUBLE, "resolutionFreqFourier2");
        MDL::addLabel(MDL_RESOLUTION_FREQREAL, LABEL_DOUBLE, "resolutionFreqReal");
        MDL::addLabel(MDL_RESOLUTION_LOCAL_RESIDUE, LABEL_DOUBLE, "localresolutionResidue");
        MDL::addLabel(MDL_RESOLUTION_LOG_STRUCTURE_FACTOR, LABEL_DOUBLE, "resolutionLogStructure");
        MDL::addLabel(MDL_RESOLUTION_STRUCTURE_FACTOR, LABEL_DOUBLE, "resolutionStructure");
        MDL::addLabel(MDL_RESOLUTION_SSNR, LABEL_DOUBLE, "resolutionSSNR");
        MDL::addLabel(MDL_RESOLUTION_RFACTOR, LABEL_DOUBLE, "resolutionRfactor");

        MDL::addLabelAlias(MDL_RESOLUTION_DPR, "DPR");
        MDL::addLabelAlias(MDL_RESOLUTION_ERRORL2, "Error_l2");
        MDL::addLabelAlias(MDL_RESOLUTION_FRC, "FRC");
        MDL::addLabelAlias(MDL_RESOLUTION_FRCRANDOMNOISE, "FRC_random_noise");
        MDL::addLabelAlias(MDL_RESOLUTION_FREQ, "Resol_Inverse_Ang");
        MDL::addLabelAlias(MDL_RESOLUTION_FREQREAL, "Resol_Ang");

        MDL::addLabel(MDL_SAMPLINGRATE, LABEL_DOUBLE, "samplingRate");
        MDL::addLabel(MDL_SAMPLINGRATE_ORIGINAL, LABEL_DOUBLE, "samplingRateOriginal");
        MDL::addLabel(MDL_SAMPLINGRATE_X, LABEL_DOUBLE, "samplingRateX");
        MDL::addLabel(MDL_SAMPLINGRATE_Y, LABEL_DOUBLE, "samplingRateY");
        MDL::addLabel(MDL_SAMPLINGRATE_Z, LABEL_DOUBLE, "samplingRateZ");

        MDL::addLabelAlias(MDL_SAMPLINGRATE, "sampling_rate"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_ORIGINAL, "sampling_rate_original"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_X, "sampling_rateX"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_Y, "sampling_rateY"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_Z, "sampling_rateZ"); //3.0

        MDL::addLabel(MDL_SCALE, LABEL_DOUBLE, "scale");
        MDL::addLabel(MDL_SCORE_BY_PCA_RESIDUAL_PROJ, LABEL_DOUBLE, "scoreByPcaResidualProj");
        MDL::addLabel(MDL_SCORE_BY_PCA_RESIDUAL_EXP, LABEL_DOUBLE, "scoreByPcaResidualExp");
        MDL::addLabel(MDL_SCORE_BY_PCA_RESIDUAL, LABEL_DOUBLE, "scoreByPcaResidual");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY, LABEL_DOUBLE, "scoreByAlignability");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_PRECISION, LABEL_DOUBLE, "scoreByAlignabilityPrecision");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_ACCURACY, LABEL_DOUBLE, "scoreByAlignabilityAccuracy");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP, LABEL_DOUBLE, "scoreByAlignabilityPrecisionExp");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF, LABEL_DOUBLE, "scoreByAlignabilityPrecisionRef");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP, LABEL_DOUBLE, "scoreByAlignabilityAccuracyExp");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF, LABEL_DOUBLE, "scoreByAlignabilityAccuracyRef");

        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_NOISE, LABEL_DOUBLE, "scoreByAlignabilityNoise");
        MDL::addLabel(MDL_SCORE_BY_EMPTINESS, LABEL_DOUBLE, "scoreEmptiness");
        MDL::addLabel(MDL_SCORE_BY_ENTROPY, LABEL_VECTOR_DOUBLE, "entropyFeatures");
        MDL::addLabel(MDL_SCORE_BY_GRANULO, LABEL_VECTOR_DOUBLE, "granuloFeatures");
        MDL::addLabel(MDL_SCORE_BY_HISTDIST, LABEL_VECTOR_DOUBLE, "histdistFeatures");
        MDL::addLabel(MDL_SCORE_BY_LBP, LABEL_VECTOR_DOUBLE, "lbpFeatures");
        MDL::addLabel(MDL_SCORE_BY_MIRROR, LABEL_DOUBLE, "scoreByMirror");
        MDL::addLabel(MDL_SCORE_BY_RAMP, LABEL_VECTOR_DOUBLE, "rampCoefficients");
        MDL::addLabel(MDL_SCORE_BY_SCREENING, LABEL_VECTOR_DOUBLE, "screenFeatures");
        MDL::addLabel(MDL_SCORE_BY_VARIANCE, LABEL_VECTOR_DOUBLE, "varianceFeatures");
        MDL::addLabel(MDL_SCORE_BY_VAR, LABEL_DOUBLE, "scoreByVariance");
        MDL::addLabel(MDL_SCORE_BY_GINI, LABEL_DOUBLE, "scoreByGiniCoeff");
        MDL::addLabel(MDL_SCORE_BY_ZERNIKE, LABEL_VECTOR_DOUBLE, "zernikeMoments");
        MDL::addLabel(MDL_SCORE_BY_ZSCORE, LABEL_DOUBLE, "scoreByZScore");

        MDL::addLabelAlias(MDL_SCALE, "Scale");
        MDL::addLabel(MDL_SELFILE, LABEL_STRING, "selfile", TAGLABEL_METADATA);
        MDL::addLabel(MDL_SERIE, LABEL_STRING, "serie");

        MDL::addLabel(MDL_SHIFT_X, LABEL_DOUBLE, "shiftX");
        MDL::addLabelAlias(MDL_SHIFT_X, "Xoff");
        MDL::addLabel(MDL_SHIFT_X2, LABEL_DOUBLE, "shiftX2");
         MDL::addLabel(MDL_SHIFT_X3, LABEL_DOUBLE, "shiftX3");
        MDL::addLabel(MDL_SHIFT_X_DIFF, LABEL_DOUBLE, "shiftXDiff");
        MDL::addLabel(MDL_SHIFT_Y, LABEL_DOUBLE, "shiftY");
        MDL::addLabelAlias(MDL_SHIFT_Y, "Yoff");
        MDL::addLabel(MDL_SHIFT_Y2, LABEL_DOUBLE, "shiftY2");
        MDL::addLabel(MDL_SHIFT_Y3, LABEL_DOUBLE, "shiftY3");
        MDL::addLabel(MDL_SHIFT_Y_DIFF, LABEL_DOUBLE, "shiftYDiff");
        MDL::addLabel(MDL_SHIFT_Z, LABEL_DOUBLE, "shiftZ");
        MDL::addLabel(MDL_SHIFT_Z2, LABEL_DOUBLE, "shiftZ2");
        MDL::addLabel(MDL_SHIFT_Z3, LABEL_DOUBLE, "shiftZ3");
        MDL::addLabelAlias(MDL_SHIFT_Z, "Zoff");
        MDL::addLabel(MDL_SHIFT_DIFF0, LABEL_DOUBLE, "shiftDiff0");
        MDL::addLabel(MDL_SHIFT_DIFF, LABEL_DOUBLE, "shiftDiff");
        MDL::addLabel(MDL_SHIFT_DIFF2, LABEL_DOUBLE, "shiftDiff2");
        MDL::addLabel(MDL_SIGMANOISE, LABEL_DOUBLE, "sigmaNoise");
        MDL::addLabel(MDL_SIGMAOFFSET, LABEL_DOUBLE, "sigmaOffset");
        MDL::addLabel(MDL_SIGNALCHANGE, LABEL_DOUBLE, "signalChange");
        MDL::addLabel(MDL_SPH_COEFFICIENTS, LABEL_VECTOR_DOUBLE, "sphCoefficients");
        MDL::addLabel(MDL_SPH_DEFORMATION, LABEL_DOUBLE, "sphDeformation");
        MDL::addLabel(MDL_SPH_TSNE_COEFF1D, LABEL_DOUBLE, "sphTsne1D");
        MDL::addLabel(MDL_SPH_TSNE_COEFF2D, LABEL_VECTOR_DOUBLE, "sphTsne2D");
        MDL::addLabel(MDL_STDDEV, LABEL_DOUBLE, "stddev");
        MDL::addLabel(MDL_STAR_COMMENT, LABEL_STRING, "starComment");
        MDL::addLabel(MDL_SUBTOMOID, LABEL_SIZET, "subtomogramId");
        MDL::addLabel(MDL_SUBTRACTION_R2, LABEL_DOUBLE, "R2subtraction");
        MDL::addLabel(MDL_SUBTRACTION_B, LABEL_DOUBLE, "Bsubtraction");
        MDL::addLabel(MDL_SUBTRACTION_BETA0, LABEL_DOUBLE, "B0subtraction");
        MDL::addLabel(MDL_SUBTRACTION_BETA1, LABEL_DOUBLE, "B1subtraction");
        MDL::addLabel(MDL_SUM, LABEL_DOUBLE, "sum");
        MDL::addLabel(MDL_SUMWEIGHT, LABEL_DOUBLE, "sumWeight");
        MDL::addLabel(MDL_SYMNO, LABEL_INT, "symNo");


		MDL::addLabel(MDL_TILTPARTICLEID, LABEL_SIZET, "tiltParticleId");
		MDL::addLabel(MDL_TOMOGRAM_VOLUME, LABEL_STRING, "tomogramVolume", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_TOMOGRAMMD, LABEL_STRING, "tomogramMetadata", TAGLABEL_METADATA);
        MDL::addLabel(MDL_TSID, LABEL_STRING, "tiltSeriesId");

        MDL::addLabel(MDL_TRANSFORM_MATRIX, LABEL_STRING, "transformMatrix");

        MDL::addLabel(MDL_TEST_SIZE, LABEL_INT, "testSize");

        MDL::addLabel(MDL_VOLUME_SCORE_SUM, LABEL_DOUBLE, "volScoreSum");
        MDL::addLabel(MDL_VOLUME_SCORE_MEAN, LABEL_DOUBLE, "volScoreMean");
        MDL::addLabel(MDL_VOLUME_SCORE_MIN, LABEL_DOUBLE, "volScoreMin");
        MDL::addLabel(MDL_VOLUME_SCORE_SUM_TH, LABEL_DOUBLE, "volScoreSumTh");
        MDL::addLabel(MDL_VOLUME_SCORE1, LABEL_DOUBLE, "volScore1");
        MDL::addLabel(MDL_VOLUME_SCORE2, LABEL_DOUBLE, "volScore2");
        MDL::addLabel(MDL_VOLUME_SCORE3, LABEL_DOUBLE, "volScore3");
        MDL::addLabel(MDL_VOLUME_SCORE4, LABEL_DOUBLE, "volScore4");
        MDL::addLabel(MDL_WEIGHT, LABEL_DOUBLE, "weight");
        MDL::addLabel(MDL_WEIGHT_P, LABEL_DOUBLE, "weight_clusterability");
        MDL::addLabelAlias(MDL_WEIGHT, "Weight");
        MDL::addLabel(MDL_WEIGHT_CONTINUOUS2, LABEL_DOUBLE, "weightContinuous2");
        MDL::addLabel(MDL_WEIGHT_JUMPER0, LABEL_DOUBLE, "weightJumper0");
        MDL::addLabel(MDL_WEIGHT_JUMPER, LABEL_DOUBLE, "weightJumper");
        MDL::addLabel(MDL_WEIGHT_JUMPER2, LABEL_DOUBLE, "weightJumper2");
        MDL::addLabel(MDL_WEIGHT_REALCORR, LABEL_DOUBLE, "weightRealCorr");
        MDL::addLabel(MDL_WEIGHT_PHASECORR, LABEL_DOUBLE, "weightPhaseCorr");
        MDL::addLabel(MDL_WEIGHT_SIGNIFICANT, LABEL_DOUBLE, "weightSignificant");
        MDL::addLabel(MDL_WEIGHT_SSNR, LABEL_DOUBLE, "weightSSNR");

        MDL::addLabel(MDL_WEIGHT_PRECISION_ALIGNABILITY, LABEL_DOUBLE, "weightPrecisionAlignability");
        MDL::addLabel(MDL_WEIGHT_ACCURACY_ALIGNABILITY, LABEL_DOUBLE, "weightAccuracyAlignability");
        MDL::addLabel(MDL_WEIGHT_ALIGNABILITY, LABEL_DOUBLE, "weightAlignability");
        MDL::addLabel(MDL_WEIGHT_PRECISION_MIRROR, LABEL_DOUBLE, "weightPrecisionMirror");

        MDL::addLabel(MDL_WROBUST, LABEL_DOUBLE, "wRobust");
        MDL::addLabel(MDL_XCOOR, LABEL_INT, "xcoor");
        MDL::addLabel(MDL_XCOOR_TILT, LABEL_INT, "xcoorTilt");

        MDL::addLabel(MDL_X, LABEL_DOUBLE, "x");
        MDL::addLabel(MDL_XSIZE, LABEL_SIZET, "xSize");
        MDL::addLabel(MDL_YCOOR, LABEL_INT, "ycoor");
        MDL::addLabel(MDL_YCOOR_TILT, LABEL_INT, "ycoorTilt");
        MDL::addLabel(MDL_Y, LABEL_DOUBLE, "y");
        MDL::addLabel(MDL_YSIZE, LABEL_SIZET, "ySize");
        MDL::addLabel(MDL_ZCOOR, LABEL_INT, "zcoor");
        MDL::addLabel(MDL_Z, LABEL_DOUBLE, "z");
        MDL::addLabel(MDL_ZSCORE, LABEL_DOUBLE, "zScore");
        MDL::addLabel(MDL_ZSCORE_HISTOGRAM, LABEL_DOUBLE, "zScoreHistogram");
        MDL::addLabel(MDL_ZSCORE_RESMEAN, LABEL_DOUBLE, "zScoreResMean");
        MDL::addLabel(MDL_ZSCORE_RESVAR, LABEL_DOUBLE, "zScoreResVar");
        MDL::addLabel(MDL_ZSCORE_RESCOV, LABEL_DOUBLE, "zScoreResCov");
        MDL::addLabel(MDL_ZSCORE_SHAPE1, LABEL_DOUBLE, "zScoreShape1");
        MDL::addLabel(MDL_ZSCORE_SHAPE2, LABEL_DOUBLE, "zScoreShape2");
        MDL::addLabel(MDL_ZSCORE_SNR1, LABEL_DOUBLE, "zScoreSNR1");
        MDL::addLabel(MDL_ZSCORE_SNR2, LABEL_DOUBLE, "zScoreSNR2");
        MDL::addLabel(MDL_ZSCORE_DEEPLEARNING1, LABEL_DOUBLE, "zScoreDeepLearning1");
        MDL::addLabel(MDL_GOOD_REGION_SCORE, LABEL_DOUBLE, "goodRegionScore");
        MDL::addLabel(MDL_ZSIZE, LABEL_SIZET, "zSize");

        MDL::addLabelAlias(MDL_XCOOR, "Xcoor");//3.0
        MDL::addLabelAlias(MDL_XCOOR, "<X position>");
        MDL::addLabelAlias(MDL_XCOOR_TILT, "XcoorTilt");//3.0

        MDL::addLabelAlias(MDL_X, "X"); //3.0
        MDL::addLabelAlias(MDL_XSIZE, "Xsize"); //3.0
        MDL::addLabelAlias(MDL_YCOOR, "Ycoor"); //3.0
        MDL::addLabelAlias(MDL_YCOOR, "<Y position>");
        MDL::addLabelAlias(MDL_YCOOR_TILT, "YcoorTilt"); //3.0
        MDL::addLabelAlias(MDL_Y, "Y"); //3.0
        MDL::addLabelAlias(MDL_YSIZE, "Ysize"); //3.0
        MDL::addLabelAlias(MDL_ZCOOR, "Zcoor"); //3.0
        MDL::addLabelAlias(MDL_Z, "Z"); //3.0
        MDL::addLabelAlias(MDL_ZSCORE, "Zscore"); //3.0
        MDL::addLabelAlias(MDL_ZSIZE, "Zsize"); //3.0


        /*Relion labels */
        MDL::addLabel(RLN_AREA_ID, LABEL_SIZET, "rlnAreaId");
        MDL::addLabel(RLN_AREA_NAME, LABEL_STRING, "rlnAreaName");

        MDL::addLabel(RLN_CTF_BFACTOR, LABEL_DOUBLE, "rlnBfactor");
        MDL::addLabel(RLN_CTF_SCALEFACTOR, LABEL_DOUBLE, "rlnCtfScalefactor");
        MDL::addLabel(RLN_CTF_VOLTAGE, LABEL_DOUBLE, "rlnVoltage");
        MDL::addLabel(RLN_CTF_DEFOCUSU, LABEL_DOUBLE, "rlnDefocusU");
        MDL::addLabel(RLN_CTF_DEFOCUSV, LABEL_DOUBLE, "rlnDefocusV");
        MDL::addLabel(RLN_CTF_DEFOCUS_ANGLE, LABEL_DOUBLE, "rlnDefocusAngle");
        MDL::addLabel(RLN_CTF_CS, LABEL_DOUBLE, "rlnSphericalAberration");
        MDL::addLabel(RLN_CTF_CA, LABEL_DOUBLE, "rlnChromaticAberration");
        MDL::addLabel(RLN_CTF_DETECTOR_PIXEL_SIZE, LABEL_DOUBLE, "rlnDetectorPixelSize");
        MDL::addLabel(RLN_CTF_ENERGY_LOSS, LABEL_DOUBLE, "rlnEnergyLoss");
        MDL::addLabel(RLN_CTF_FOM, LABEL_DOUBLE, "rlnCtfFigureOfMerit");
        MDL::addLabel(RLN_CTF_IMAGE, LABEL_STRING, "rlnCtfImage", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_CTF_LENS_STABILITY, LABEL_DOUBLE, "rlnLensStability");
        MDL::addLabel(RLN_CTF_MAGNIFICATION, LABEL_DOUBLE, "rlnMagnification");
        MDL::addLabel(RLN_CTF_MAXRES, LABEL_DOUBLE, "rlnCtfMaxResolution");
        MDL::addLabel(RLN_CTF_CONVERGENCE_CONE, LABEL_DOUBLE, "rlnConvergenceCone");
        MDL::addLabel(RLN_CTF_LONGITUDINAL_DISPLACEMENT, LABEL_DOUBLE, "rlnLongitudinalDisplacement");
        MDL::addLabel(RLN_CTF_TRANSVERSAL_DISPLACEMENT, LABEL_DOUBLE, "rlnTransversalDisplacement");
        MDL::addLabel(RLN_CTF_Q0, LABEL_DOUBLE, "rlnAmplitudeContrast");
        MDL::addLabel(RLN_CTF_VALIDATIONSCORE, LABEL_DOUBLE, "rlnCtfValidationScore");
        MDL::addLabel(RLN_CTF_VALUE, LABEL_DOUBLE, "rlnCtfValue");
        MDL::addLabel(RLN_CTF_PHASESHIFT, LABEL_DOUBLE, "rlnPhaseShift");

        MDL::addLabel(RLN_IMAGE_NAME, LABEL_STRING, "rlnImageName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_RECONSTRUCT_NAME, LABEL_STRING, "rlnReconstructImageName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_ID, LABEL_SIZET, "rlnImageId");
        MDL::addLabel(RLN_IMAGE_ENABLED, LABEL_BOOL, "rlnEnabled");
        MDL::addLabel(RLN_IMAGE_DATATYPE, LABEL_INT, "rlnDataType");
        // Label rlnImageDimensionality is originally rlnDataDimensionality, which is
        // duplicated for other label. A relion bug???
        MDL::addLabel(RLN_IMAGE_DIMENSIONALITY, LABEL_INT, "rlnImageDimensionality");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_X, LABEL_DOUBLE, "rlnBeamTiltX");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_Y, LABEL_DOUBLE, "rlnBeamTiltY");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_GROUP, LABEL_STRING, "rlnBeamTiltGroupName");
        MDL::addLabel(RLN_IMAGE_COORD_X, LABEL_DOUBLE, "rlnCoordinateX");
        MDL::addLabel(RLN_IMAGE_COORD_Y, LABEL_DOUBLE, "rlnCoordinateY");
        MDL::addLabel(RLN_IMAGE_COORD_Z, LABEL_DOUBLE, "rlnCoordinateZ");
        MDL::addLabel(RLN_IMAGE_FRAME_NR, LABEL_INT, "rlnMovieFrameNumber");
        MDL::addLabel(RLN_IMAGE_MAGNIFICATION_CORRECTION, LABEL_DOUBLE, "rlnMagnificationCorrection");
        MDL::addLabel(RLN_IMAGE_NORM_CORRECTION, LABEL_DOUBLE, "rlnNormCorrection");
        MDL::addLabel(RLN_IMAGE_ORI_NAME, LABEL_STRING, "rlnImageOriginalName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE, LABEL_DOUBLE, "rlnSamplingRate");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_X, LABEL_DOUBLE, "rlnSamplingRateX");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_Y, LABEL_DOUBLE, "rlnSamplingRateY");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_Z, LABEL_DOUBLE, "rlnSamplingRateZ");
        MDL::addLabel(RLN_IMAGE_SIZE, LABEL_INT, "rlnImageSize");
        MDL::addLabel(RLN_IMAGE_SIZEX, LABEL_INT, "rlnImageSizeX");
        MDL::addLabel(RLN_IMAGE_SIZEY, LABEL_INT, "rlnImageSizeY");
        MDL::addLabel(RLN_IMAGE_SIZEZ, LABEL_INT, "rlnImageSizeZ");
        MDL::addLabel(RLN_IMAGE_STATS_MIN, LABEL_DOUBLE, "rlnMinimumValue");
        MDL::addLabel(RLN_IMAGE_STATS_MAX, LABEL_DOUBLE, "rlnMaximumValue");
        MDL::addLabel(RLN_IMAGE_STATS_AVG, LABEL_DOUBLE, "rlnAverageValue");
        MDL::addLabel(RLN_IMAGE_STATS_STDDEV, LABEL_DOUBLE, "rlnStandardDeviationValue");
        MDL::addLabel(RLN_IMAGE_STATS_SKEW, LABEL_DOUBLE, "rlnSkewnessValue");
        MDL::addLabel(RLN_IMAGE_STATS_KURT, LABEL_DOUBLE, "rlnKurtosisExcessValue");
        MDL::addLabel(RLN_IMAGE_WEIGHT, LABEL_DOUBLE, "rlnImageWeight");

        MDL::addLabel(RLN_MASK_NAME, LABEL_STRING, "rlnMaskName");

        MDL::addLabel(RLN_MATRIX_1_1, LABEL_DOUBLE, "rlnMatrix_1_1");
        MDL::addLabel(RLN_MATRIX_1_2, LABEL_DOUBLE, "rlnMatrix_1_2");
        MDL::addLabel(RLN_MATRIX_1_3, LABEL_DOUBLE, "rlnMatrix_1_3");
        MDL::addLabel(RLN_MATRIX_2_1, LABEL_DOUBLE, "rlnMatrix_2_1");
        MDL::addLabel(RLN_MATRIX_2_2, LABEL_DOUBLE, "rlnMatrix_2_2");
        MDL::addLabel(RLN_MATRIX_2_3, LABEL_DOUBLE, "rlnMatrix_2_3");
        MDL::addLabel(RLN_MATRIX_3_1, LABEL_DOUBLE, "rlnMatrix_3_1");
        MDL::addLabel(RLN_MATRIX_3_2, LABEL_DOUBLE, "rlnMatrix_3_2");
        MDL::addLabel(RLN_MATRIX_3_3, LABEL_DOUBLE, "rlnMatrix_3_3");

        MDL::addLabel(RLN_MICROGRAPH_ID, LABEL_SIZET, "rlnMicrographId");
        MDL::addLabel(RLN_MICROGRAPH_MOVIE_NAME, LABEL_STRING, "rlnMicrographMovieName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MICROGRAPH_NAME, LABEL_STRING, "rlnMicrographName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MICROGRAPH_NAME_WODOSE, LABEL_STRING, "rlnMicrographNameNoDW", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MICROGRAPH_TILT_ANGLE, LABEL_DOUBLE, "rlnMicrographTiltAngle");
        MDL::addLabel(RLN_MICROGRAPH_TILT_AXIS_DIRECTION, LABEL_DOUBLE, "rlnMicrographTiltAxisDirection");
        MDL::addLabel(RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE, LABEL_DOUBLE, "rlnMicrographTiltAxisOutOfPlane");

        MDL::addLabel(RLN_MLMODEL_ACCURACY_ROT, LABEL_DOUBLE, "rlnAccuracyRotations");
        MDL::addLabel(RLN_MLMODEL_ACCURACY_TRANS, LABEL_DOUBLE, "rlnAccuracyTranslations");
        MDL::addLabel(RLN_MLMODEL_AVE_PMAX, LABEL_DOUBLE, "rlnAveragePmax");
        MDL::addLabel(RLN_MLMODEL_CURRENT_RESOLUTION, LABEL_DOUBLE, "rlnCurrentResolution");
        MDL::addLabel(RLN_MLMODEL_CURRENT_SIZE, LABEL_INT, "rlnCurrentImageSize");
        MDL::addLabel(RLN_MLMODEL_DATA_VS_PRIOR_REF, LABEL_DOUBLE, "rlnSsnrMap");
        MDL::addLabel(RLN_MLMODEL_DIMENSIONALITY, LABEL_INT, "rlnReferenceDimensionality");
        MDL::addLabel(RLN_MLMODEL_DIMENSIONALITY_DATA, LABEL_INT, "rlnDataDimensionality");
        MDL::addLabel(RLN_MLMODEL_DIFF2_HALVES_REF, LABEL_DOUBLE, "rlnDiff2RandomHalves");
        MDL::addLabel(RLN_MLMODEL_ESTIM_RESOL_REF, LABEL_DOUBLE, "rlnEstimatedResolution");
        MDL::addLabel(RLN_MLMODEL_FOURIER_COVERAGE_REF, LABEL_DOUBLE, "rlnFourierCompleteness");
        MDL::addLabel(RLN_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, LABEL_DOUBLE, "rlnOverallFourierCompleteness");
        MDL::addLabel(RLN_MLMODEL_FSC_HALVES_REF, LABEL_DOUBLE, "rlnGoldStandardFsc");
        MDL::addLabel(RLN_MLMODEL_GROUP_NAME, LABEL_STRING, "rlnGroupName");
        MDL::addLabel(RLN_MLMODEL_GROUP_NO, LABEL_SIZET, "rlnGroupNumber");
        MDL::addLabel(RLN_MLMODEL_GROUP_NR_PARTICLES, LABEL_SIZET, "rlnGroupNrParticles");
        MDL::addLabel(RLN_MLMODEL_GROUP_SCALE_CORRECTION, LABEL_DOUBLE, "rlnGroupScaleCorrection");
        MDL::addLabel(RLN_MLMODEL_HELICAL_NR_ASU, LABEL_INT, "rlnNrHelicalAsymUnits");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST, LABEL_DOUBLE, "rlnHelicalTwist");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST_MIN, LABEL_DOUBLE, "rlnHelicalTwistMin");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST_MAX, LABEL_DOUBLE, "rlnHelicalTwistMax");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST_INITIAL_STEP, LABEL_DOUBLE, "rlnHelicalTwistInitialStep");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE, LABEL_DOUBLE, "rlnHelicalRise");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE_MIN, LABEL_DOUBLE, "rlnHelicalRiseMin");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE_MAX, LABEL_DOUBLE, "rlnHelicalRiseMax");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE_INITIAL_STEP, LABEL_DOUBLE, "rlnHelicalRiseInitialStep");
        MDL::addLabel(RLN_MLMODEL_INTERPOLATOR, LABEL_INT, "rlnFourierSpaceInterpolator");
        MDL::addLabel(RLN_MLMODEL_IS_HELIX, LABEL_BOOL, "rlnIsHelix");
        MDL::addLabel(RLN_MLMODEL_LL, LABEL_DOUBLE, "rlnLogLikelihood");
        MDL::addLabel(RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,  LABEL_INT, "rlnMinRadiusNnInterpolation");
        MDL::addLabel(RLN_MLMODEL_NORM_CORRECTION_AVG, LABEL_DOUBLE, "rlnNormCorrectionAverage");
        MDL::addLabel(RLN_MLMODEL_NR_BODIES, LABEL_INT, "rlnNrBodies");
        MDL::addLabel(RLN_MLMODEL_NR_CLASSES, LABEL_INT, "rlnNrClasses");
        MDL::addLabel(RLN_MLMODEL_NR_GROUPS, LABEL_INT, "rlnNrGroups");
        MDL::addLabel(RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION, LABEL_DOUBLE, "rlnSpectralOrientabilityContribution");
        MDL::addLabel(RLN_MLMODEL_ORIGINAL_SIZE, LABEL_INT, "rlnOriginalImageSize");
        MDL::addLabel(RLN_MLMODEL_PADDING_FACTOR, LABEL_DOUBLE, "rlnPaddingFactor");
        MDL::addLabel(RLN_MLMODEL_PDF_CLASS, LABEL_DOUBLE, "rlnClassDistribution");
        MDL::addLabel(RLN_MLMODEL_PRIOR_OFFX_CLASS, LABEL_DOUBLE, "rlnClassPriorOffsetX");
        MDL::addLabel(RLN_MLMODEL_PRIOR_OFFY_CLASS, LABEL_DOUBLE, "rlnClassPriorOffsetY");
        MDL::addLabel(RLN_MLMODEL_PDF_ORIENT, LABEL_DOUBLE, "rlnOrientationDistribution");
        MDL::addLabel(RLN_MLMODEL_PIXEL_SIZE, LABEL_DOUBLE, "rlnPixelSize");
        MDL::addLabel(RLN_MLMODEL_POWER_REF, LABEL_DOUBLE, "rlnReferenceSpectralPower");
        MDL::addLabel(RLN_MLMODEL_PRIOR_MODE, LABEL_INT, "rlnOrientationalPriorMode");
        MDL::addLabel(RLN_MLMODEL_REF_IMAGE, LABEL_STRING, "rlnReferenceImage", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MLMODEL_SGD_GRADIENT_IMAGE, LABEL_STRING, "rlnSGDGradientImage");
        MDL::addLabel(RLN_MLMODEL_SIGMA_OFFSET, LABEL_DOUBLE, "rlnSigmaOffsets");
        MDL::addLabel(RLN_MLMODEL_SIGMA2_NOISE, LABEL_DOUBLE, "rlnSigma2Noise");
        MDL::addLabel(RLN_MLMODEL_SIGMA2_REF, LABEL_DOUBLE, "rlnReferenceSigma2");
        MDL::addLabel(RLN_MLMODEL_SIGMA_ROT, LABEL_DOUBLE, "rlnSigmaPriorRotAngle");
        MDL::addLabel(RLN_MLMODEL_SIGMA_TILT, LABEL_DOUBLE, "rlnSigmaPriorTiltAngle");
        MDL::addLabel(RLN_MLMODEL_SIGMA_PSI, LABEL_DOUBLE, "rlnSigmaPriorPsiAngle");
        MDL::addLabel(RLN_MLMODEL_SSNR_REF, LABEL_DOUBLE, "rlnSignalToNoiseRatio");
        MDL::addLabel(RLN_MLMODEL_TAU2_FUDGE_FACTOR, LABEL_DOUBLE, "rlnTau2FudgeFactor");
        MDL::addLabel(RLN_MLMODEL_TAU2_REF, LABEL_DOUBLE, "rlnReferenceTau2");

        MDL::addLabel(RLN_OPTIMISER_ACCURACY_ROT, LABEL_DOUBLE, "rlnOverallAccuracyRotations");
        MDL::addLabel(RLN_OPTIMISER_ACCURACY_TRANS, LABEL_DOUBLE, "rlnOverallAccuracyTranslations");
        MDL::addLabel(RLN_OPTIMISER_ADAPTIVE_FRACTION, LABEL_DOUBLE, "rlnAdaptiveOversampleFraction");
        MDL::addLabel(RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING, LABEL_INT, "rlnAdaptiveOversampleOrder");
        MDL::addLabel(RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER, LABEL_INT, "rlnAutoLocalSearchesHealpixOrder");
        MDL::addLabel(RLN_OPTIMISER_AVAILABLE_MEMORY, LABEL_DOUBLE, "rlnAvailableMemory");
        MDL::addLabel(RLN_OPTIMISER_BEST_RESOL_THUS_FAR, LABEL_DOUBLE, "rlnBestResolutionThusFar");
        MDL::addLabel(RLN_OPTIMISER_COARSE_SIZE, LABEL_INT, "rlnCoarseImageSize");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, LABEL_DOUBLE, "rlnChangesOptimalOffsets");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, LABEL_DOUBLE, "rlnChangesOptimalOrientations");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES, LABEL_DOUBLE, "rlnChangesOptimalClasses");
        MDL::addLabel(RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, LABEL_BOOL, "rlnCtfDataArePhaseFlipped");
        MDL::addLabel(RLN_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, LABEL_BOOL, "rlnCtfDataAreCtfPremultiplied");
        MDL::addLabel(RLN_OPTIMISER_DATA_STARFILE, LABEL_STRING, "rlnExperimentalDataStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_DO_AUTO_REFINE, LABEL_BOOL, "rlnDoAutoRefine");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_CTF, LABEL_BOOL, "rlnDoCorrectCtf");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION, LABEL_BOOL, "rlnDoCorrectMagnification");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_NORM, LABEL_BOOL, "rlnDoCorrectNorm");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_SCALE, LABEL_BOOL, "rlnDoCorrectScale");
        MDL::addLabel(RLN_OPTIMISER_DO_HELICAL_REFINE, LABEL_BOOL, "rlnDoHelicalRefine");
        MDL::addLabel(RLN_OPTIMISER_DO_MAP, LABEL_BOOL, "rlnDoMapEstimation");
        MDL::addLabel(RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, LABEL_BOOL, "rlnDoOnlyFlipCtfPhases");
        MDL::addLabel(RLN_OPTIMISER_DO_REALIGN_MOVIES, LABEL_BOOL, "rlnDoRealignMovies");
        MDL::addLabel(RLN_OPTIMISER_DO_SGD, LABEL_BOOL, "rlnDoStochasticGradientDescent");
        MDL::addLabel(RLN_OPTIMISER_DO_SOLVENT_FLATTEN, LABEL_BOOL, "rlnDoSolventFlattening");
        MDL::addLabel(RLN_OPTIMISER_DO_SKIP_ALIGN, LABEL_BOOL, "rlnDoSkipAlign");
        MDL::addLabel(RLN_OPTIMISER_DO_SKIP_ROTATE, LABEL_BOOL, "rlnDoSkipRotate");
        MDL::addLabel(RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES, LABEL_BOOL, "rlnDoSplitRandomHalves");
        MDL::addLabel(RLN_OPTIMISER_DO_ZERO_MASK, LABEL_BOOL, "rlnDoZeroMask");
        MDL::addLabel(RLN_OPTIMISER_FIX_SIGMA_NOISE, LABEL_BOOL, "rlnFixSigmaNoiseEstimates");
        MDL::addLabel(RLN_OPTIMISER_FIX_SIGMA_OFFSET ,LABEL_BOOL, "rlnFixSigmaOffsetEstimates");
        MDL::addLabel(RLN_OPTIMISER_FIX_TAU, LABEL_BOOL, "rlnFixTauEstimates");
        MDL::addLabel(RLN_OPTIMISER_HAS_CONVERGED, LABEL_BOOL, "rlnHasConverged");
        MDL::addLabel(RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, LABEL_BOOL, "rlnHasHighFscAtResolLimit");
        MDL::addLabel(RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, LABEL_INT, "rlnHasLargeSizeIncreaseIterationsAgo");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_TWIST_INITIAL, LABEL_DOUBLE, "rlnHelicalTwistInitial");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_RISE_INITIAL, LABEL_DOUBLE, "rlnHelicalRiseInitial");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_Z_PERCENTAGE, LABEL_DOUBLE, "rlnHelicalCentralProportion");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, LABEL_DOUBLE, "rlnHelicalMaskTubeInnerDiameter");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, LABEL_DOUBLE, "rlnHelicalMaskTubeOuterDiameter");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, LABEL_BOOL, "rlnHelicalSymmetryLocalRefinement");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_SIGMA_DISTANCE, LABEL_DOUBLE, "rlnHelicalSigmaDistance");
        MDL::addLabel(RLN_OPTIMISER_IGNORE_HELICAL_SYMMETRY, LABEL_BOOL, "rlnIgnoreHelicalSymmetry");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, LABEL_BOOL, "rlnHelicalKeepTiltPriorFixed");
        MDL::addLabel(RLN_OPTIMISER_HIGHRES_LIMIT_EXP, LABEL_DOUBLE, "rlnHighresLimitExpectation");
        MDL::addLabel(RLN_OPTIMISER_HIGHRES_LIMIT_SGD, LABEL_DOUBLE, "rlnHighresLimitSGD");
        MDL::addLabel(RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, LABEL_BOOL, "rlnDoIgnoreCtfUntilFirstPeak");
        MDL::addLabel(RLN_OPTIMISER_INCR_SIZE, LABEL_INT, "rlnIncrementImageSize");
        MDL::addLabel(RLN_OPTIMISER_ITERATION_NO, LABEL_INT, "rlnCurrentIteration");
        MDL::addLabel(RLN_OPTIMISER_LOCAL_SYMMETRY_FILENAME, LABEL_STRING, "rlnLocalSymmetryFile");
        MDL::addLabel(RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, LABEL_DOUBLE, "rlnJoinHalvesUntilThisResolution");
        MDL::addLabel(RLN_OPTIMISER_MAGNIFICATION_RANGE, LABEL_DOUBLE, "rlnMagnificationSearchRange");
        MDL::addLabel(RLN_OPTIMISER_MAGNIFICATION_STEP, LABEL_DOUBLE, "rlnMagnificationSearchStep");
        MDL::addLabel(RLN_OPTIMISER_MAX_COARSE_SIZE, LABEL_INT, "rlnMaximumCoarseImageSize");
        MDL::addLabel(RLN_OPTIMISER_MAX_NR_POOL, LABEL_INT, "rlnMaxNumberOfPooledParticles");
        MDL::addLabel(RLN_OPTIMISER_MODEL_STARFILE, LABEL_STRING, "rlnModelStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_MODEL_STARFILE2, LABEL_STRING, "rlnModelStarFile2", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_NR_ITERATIONS, LABEL_INT, "rlnNumberOfIterations");
        MDL::addLabel(RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN, LABEL_INT, "rlnNumberOfIterWithoutResolutionGain");
        MDL::addLabel(RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, LABEL_INT, "rlnNumberOfIterWithoutChangingAssignments");
        MDL::addLabel(RLN_OPTIMISER_OUTPUT_ROOTNAME, LABEL_STRING, "rlnOutputRootName");
        MDL::addLabel(RLN_OPTIMISER_PARTICLE_DIAMETER, LABEL_DOUBLE, "rlnParticleDiameter");
        MDL::addLabel(RLN_OPTIMISER_RADIUS_MASK_3D_MAP, LABEL_INT, "rlnRadiusMaskMap");
        MDL::addLabel(RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES, LABEL_INT, "rlnRadiusMaskExpImages");
        MDL::addLabel(RLN_OPTIMISER_RANDOM_SEED, LABEL_INT, "rlnRandomSeed");
        MDL::addLabel(RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED, LABEL_BOOL, "rlnRefsAreCtfCorrected");
        MDL::addLabel(RLN_OPTIMISER_SGD_MU, LABEL_DOUBLE, "rlnSgdMuFactor");
        MDL::addLabel(RLN_OPTIMISER_SGD_SIGMA2FUDGE_INI, LABEL_DOUBLE, "rlnSgdSigma2FudgeInitial");
        MDL::addLabel(RLN_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, LABEL_SIZET, "rlnSgdSigma2FudgeHalflife");
        MDL::addLabel(RLN_OPTIMISER_SGD_SUBSET_START, LABEL_INT, "rlnSgdNextSubset");
        MDL::addLabel(RLN_OPTIMISER_SGD_SUBSET_SIZE, LABEL_SIZET, "rlnSgdSubsetSize");
        MDL::addLabel(RLN_OPTIMISER_SGD_WRITE_EVERY_SUBSET, LABEL_INT, "rlnSgdWriteEverySubset");
        MDL::addLabel(RLN_OPTIMISER_SGD_MAX_SUBSETS, LABEL_SIZET, "rlnSgdMaxSubsets");
        MDL::addLabel(RLN_OPTIMISER_SGD_STEPSIZE, LABEL_DOUBLE, "rlnSgdStepsize");


        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, LABEL_INT, "rlnSmallestChangesClasses");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, LABEL_DOUBLE, "rlnSmallestChangesOffsets");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, LABEL_DOUBLE, "rlnSmallestChangesOrientations");
        MDL::addLabel(RLN_OPTIMISER_SAMPLING_STARFILE, LABEL_STRING, "rlnOrientSamplingStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_SOLVENT_MASK_NAME, LABEL_STRING, "rlnSolventMaskName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_OPTIMISER_SOLVENT_MASK2_NAME, LABEL_STRING, "rlnSolventMask2Name", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_OPTIMISER_TAU_SPECTRUM_NAME, LABEL_STRING, "rlnTauSpectrumName");
        MDL::addLabel(RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING, LABEL_BOOL, "rlnUseTooCoarseSampling");
        MDL::addLabel(RLN_OPTIMISER_WIDTH_MASK_EDGE, LABEL_INT, "rlnWidthMaskEdge");

        MDL::addLabel(RLN_ORIENT_FLIP, LABEL_BOOL, "rlnIsFlip");
        MDL::addLabel(RLN_ORIENT_ID, LABEL_SIZET, "rlnOrientationsID");
        MDL::addLabel(RLN_ORIENT_ORIGIN_X, LABEL_DOUBLE, "rlnOriginX");
        MDL::addLabel(RLN_ORIENT_ORIGIN_X_PRIOR, LABEL_DOUBLE, "rlnOriginXPrior");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Y, LABEL_DOUBLE, "rlnOriginY");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Y_PRIOR, LABEL_DOUBLE, "rlnOriginYPrior");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Z, LABEL_DOUBLE, "rlnOriginZ");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Z_PRIOR, LABEL_DOUBLE, "rlnOriginZPrior");
        MDL::addLabel(RLN_ORIENT_ROT, LABEL_DOUBLE, "rlnAngleRot");
        MDL::addLabel(RLN_ORIENT_ROT_PRIOR, LABEL_DOUBLE, "rlnAngleRotPrior");
        MDL::addLabel(RLN_ORIENT_TILT, LABEL_DOUBLE, "rlnAngleTilt");
        MDL::addLabel(RLN_ORIENT_TILT_PRIOR, LABEL_DOUBLE, "rlnAngleTiltPrior");
        MDL::addLabel(RLN_ORIENT_PSI, LABEL_DOUBLE, "rlnAnglePsi");
        MDL::addLabel(RLN_ORIENT_PSI_PRIOR, LABEL_DOUBLE, "rlnAnglePsiPrior");
        MDL::addLabel(RLN_ORIENT_PSI_PRIOR_FLIP_RATIO, LABEL_DOUBLE, "rlnAnglePsiFlipRatio");

        MDL::addLabel(RLN_PARTICLE_AUTOPICK_FOM, LABEL_DOUBLE, "rlnAutopickFigureOfMerit");
        MDL::addLabel(RLN_PARTICLE_CLASS, LABEL_INT, "rlnClassNumber");
        MDL::addLabel(RLN_PARTICLE_DLL, LABEL_DOUBLE, "rlnLogLikeliContribution");
        MDL::addLabel(RLN_PARTICLE_ID, LABEL_SIZET, "rlnParticleId");
        MDL::addLabel(RLN_PARTICLE_FOM, LABEL_DOUBLE, "rlnParticleFigureOfMerit");
        MDL::addLabel(RLN_PARTICLE_HELICAL_TUBE_ID, LABEL_INT, "rlnHelicalTubeID");
        MDL::addLabel(RLN_PARTICLE_HELICAL_TUBE_PITCH, LABEL_DOUBLE, "rlnHelicalTubePitch");
        MDL::addLabel(RLN_PARTICLE_HELICAL_TRACK_LENGTH, LABEL_DOUBLE, "rlnHelicalTrackLength");
        MDL::addLabel(RLN_PARTICLE_KL_DIVERGENCE, LABEL_DOUBLE, "rlnKullbackLeibnerDivergence");
        MDL::addLabel(RLN_PARTICLE_MOVIE_RUNNING_AVG, LABEL_INT, "rlnMovieFramesRunningAverage");
        MDL::addLabel(RLN_PARTICLE_NAME, LABEL_STRING, "rlnParticleName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES, LABEL_INT, "rlnNrOfSignificantSamples"); /**< particle, Number of orientations contributing to weights*/
        MDL::addLabel(RLN_PARTICLE_NR_FRAMES, LABEL_INT, "rlnNrOfFrames");
        MDL::addLabel(RLN_PARTICLE_NR_FRAMES_AVG, LABEL_INT, "rlnAverageNrOfFrames");
        MDL::addLabel(RLN_PARTICLE_ORI_NAME, LABEL_STRING, "rlnOriginalParticleName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_PARTICLE_PMAX, LABEL_DOUBLE, "rlnMaxValueProbDistribution"); /**< particle, Maximum value of probability distribution */
        MDL::addLabel(RLN_PARTICLE_RANDOM_SUBSET, LABEL_INT, "rlnRandomSubset");

        MDL::addLabel(RLN_PIPELINE_EDGE_FROM, LABEL_STRING , "rlnPipeLineEdgeFromNode");
        MDL::addLabel(RLN_PIPELINE_EDGE_TO, LABEL_STRING ,"rlnPipeLineEdgeToNode");
        MDL::addLabel(RLN_PIPELINE_EDGE_PROCESS, LABEL_STRING ,"rlnPipeLineEdgeProcess");
        MDL::addLabel(RLN_PIPELINE_JOB_COUNTER, LABEL_INT, "rlnPipeLineJobCounter");
        MDL::addLabel(RLN_PIPELINE_NODE_NAME, LABEL_STRING , "rlnPipeLineNodeName");
        MDL::addLabel(RLN_PIPELINE_NODE_TYPE, LABEL_INT, "rlnPipeLineNodeType");
        MDL::addLabel(RLN_PIPELINE_PROCESS_ALIAS, LABEL_STRING , "rlnPipeLineProcessAlias");
        MDL::addLabel(RLN_PIPELINE_PROCESS_NAME, LABEL_STRING , "rlnPipeLineProcessName");
        MDL::addLabel(RLN_PIPELINE_PROCESS_TYPE, LABEL_INT, "rlnPipeLineProcessType");
        MDL::addLabel(RLN_PIPELINE_PROCESS_STATUS, LABEL_INT, "rlnPipeLineProcessStatus");

        MDL::addLabel(RLN_POSTPROCESS_AMPLCORR_MASKED, LABEL_DOUBLE, "rlnAmplitudeCorrelationMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_AMPLCORR_UNMASKED,  LABEL_DOUBLE, "rlnAmplitudeCorrelationUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_BFACTOR, LABEL_DOUBLE, "rlnBfactorUsedForSharpening");
        MDL::addLabel(RLN_POSTPROCESS_DPR_MASKED, LABEL_DOUBLE, "rlnDifferentialPhaseResidualMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_DPR_UNMASKED,  LABEL_DOUBLE, "rlnDifferentialPhaseResidualUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FINAL_RESOLUTION, LABEL_DOUBLE, "rlnFinalResolution");
        MDL::addLabel(RLN_POSTPROCESS_FSC_GENERAL, LABEL_DOUBLE, "rlnFourierShellCorrelation");
        MDL::addLabel(RLN_POSTPROCESS_FSC_TRUE, LABEL_DOUBLE, "rlnFourierShellCorrelationCorrected");
        MDL::addLabel(RLN_POSTPROCESS_FSC_MASKED, LABEL_DOUBLE, "rlnFourierShellCorrelationMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FSC_UNMASKED, LABEL_DOUBLE, "rlnFourierShellCorrelationUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FSC_RANDOM_MASKED, LABEL_DOUBLE, "rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT, LABEL_DOUBLE, "rlnFittedInterceptGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_SLOPE, LABEL_DOUBLE, "rlnFittedSlopeGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_CORRELATION, LABEL_DOUBLE, "rlnCorrelationFitGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_IN, LABEL_DOUBLE, "rlnLogAmplitudesOriginal");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_INVMTF, LABEL_DOUBLE, "rlnLogAmplitudesMTFCorrected");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED, LABEL_DOUBLE, "rlnLogAmplitudesWeighted");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED, LABEL_DOUBLE, "rlnLogAmplitudesSharpened");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT, LABEL_DOUBLE, "rlnLogAmplitudesIntercept");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_RESOL_SQUARED, LABEL_DOUBLE, "rlnResolutionSquared");
        MDL::addLabel(RLN_POSTPROCESS_MTF_VALUE, LABEL_DOUBLE, "rlnMtfValue");

        MDL::addLabel(RLN_SAMPLING_IS_3D, LABEL_BOOL, "rlnIs3DSampling");
        MDL::addLabel(RLN_SAMPLING_IS_3D_TRANS, LABEL_BOOL, "rlnIs3DTranslationalSampling");
        MDL::addLabel(RLN_SAMPLING_HEALPIX_ORDER, LABEL_INT, "rlnHealpixOrder");
        MDL::addLabel(RLN_SAMPLING_HELICAL_OFFSET_STEP, LABEL_DOUBLE, "rlnHelicalOffsetStep");
        MDL::addLabel(RLN_SAMPLING_LIMIT_TILT, LABEL_DOUBLE, "rlnTiltAngleLimit");
        MDL::addLabel(RLN_SAMPLING_OFFSET_RANGE, LABEL_DOUBLE, "rlnOffsetRange");
        MDL::addLabel(RLN_SAMPLING_OFFSET_STEP, LABEL_DOUBLE, "rlnOffsetStep");
        MDL::addLabel(RLN_SAMPLING_PERTURB, LABEL_DOUBLE, "rlnSamplingPerturbInstance");
        MDL::addLabel(RLN_SAMPLING_PERTURBATION_FACTOR, LABEL_DOUBLE, "rlnSamplingPerturbFactor");
        MDL::addLabel(RLN_SAMPLING_PSI_STEP, LABEL_DOUBLE, "rlnPsiStep");
        MDL::addLabel(RLN_SAMPLING_SYMMETRY, LABEL_STRING, "rlnSymmetryGroup");

        MDL::addLabel(RLN_SELECTED, LABEL_BOOL, "rlnSelected");
        MDL::addLabel(RLN_SELECT_PARTICLES_ZSCORE, LABEL_DOUBLE, "rlnParticleSelectZScore");
        MDL::addLabel(RLN_SORTED_IDX, LABEL_SIZET, "rlnSortedIndex");
        MDL::addLabel(RLN_STARFILE_MOVIE_PARTICLES, LABEL_STRING, "rlnStarFileMovieParticles");
        MDL::addLabel(RLN_PERFRAME_CUMULATIVE_WEIGHT, LABEL_DOUBLE, "rlnPerFrameCumulativeWeight");
        MDL::addLabel(RLN_PERFRAME_RELATIVE_WEIGHT, LABEL_DOUBLE, "rlnPerFrameRelativeWeight");
        MDL::addLabel(RLN_RESOLUTION, LABEL_DOUBLE, "rlnResolution");
        MDL::addLabel(RLN_RESOLUTION_ANGSTROM, LABEL_DOUBLE, "rlnAngstromResolution");
        MDL::addLabel(RLN_RESOLUTION_INVPIXEL, LABEL_DOUBLE, "rlnResolutionInversePixel");
        MDL::addLabel(RLN_SPECTRAL_IDX, LABEL_INT, "rlnSpectralIndex");

        MDL::addLabelAlias(RLN_CTF_BFACTOR, "rlnCtfBfactor"); //Relion-2.0

        /** Buffer labels */
        MDL::addLabel(BUFFER_01, LABEL_STRING, "buffer_01");
        MDL::addLabel(BUFFER_02, LABEL_STRING, "buffer_02");
        MDL::addLabel(BUFFER_03, LABEL_STRING, "buffer_03");
        MDL::addLabel(BUFFER_04, LABEL_STRING, "buffer_04");
        MDL::addLabel(BUFFER_05, LABEL_STRING, "buffer_05");
        MDL::addLabel(BUFFER_06, LABEL_STRING, "buffer_06");
        MDL::addLabel(BUFFER_07, LABEL_STRING, "buffer_07");
        MDL::addLabel(BUFFER_08, LABEL_STRING, "buffer_08");
        MDL::addLabel(BUFFER_09, LABEL_STRING, "buffer_09");
        MDL::addLabel(BUFFER_10, LABEL_STRING, "buffer_10");
        MDL::addLabel(BUFFER_11, LABEL_STRING, "buffer_11");
        MDL::addLabel(BUFFER_12, LABEL_STRING, "buffer_12");
        MDL::addLabel(BUFFER_13, LABEL_STRING, "buffer_13");
        MDL::addLabel(BUFFER_14, LABEL_STRING, "buffer_14");
        MDL::addLabel(BUFFER_15, LABEL_STRING, "buffer_15");
        MDL::addLabel(BUFFER_16, LABEL_STRING, "buffer_16");
        MDL::addLabel(BUFFER_17, LABEL_STRING, "buffer_17");
        MDL::addLabel(BUFFER_18, LABEL_STRING, "buffer_18");
        MDL::addLabel(BUFFER_19, LABEL_STRING, "buffer_19");
        MDL::addLabel(BUFFER_20, LABEL_STRING, "buffer_20");
        MDL::addLabel(BUFFER_21, LABEL_STRING, "buffer_21");
        MDL::addLabel(BUFFER_22, LABEL_STRING, "buffer_22");
        MDL::addLabel(BUFFER_23, LABEL_STRING, "buffer_23");
        MDL::addLabel(BUFFER_24, LABEL_STRING, "buffer_24");
        MDL::addLabel(BUFFER_25, LABEL_STRING, "buffer_25");
        MDL::addLabel(BUFFER_26, LABEL_STRING, "buffer_26");
        MDL::addLabel(BUFFER_27, LABEL_STRING, "buffer_27");
        MDL::addLabel(BUFFER_28, LABEL_STRING, "buffer_28");
        MDL::addLabel(BUFFER_29, LABEL_STRING, "buffer_29");
        MDL::addLabel(BUFFER_30, LABEL_STRING, "buffer_30");
        MDL::addLabel(BUFFER_31, LABEL_STRING, "buffer_31");
        MDL::addLabel(BUFFER_32, LABEL_STRING, "buffer_32");
        MDL::addLabel(BUFFER_33, LABEL_STRING, "buffer_33");
        MDL::addLabel(BUFFER_34, LABEL_STRING, "buffer_34");
        MDL::addLabel(BUFFER_35, LABEL_STRING, "buffer_35");
        MDL::addLabel(BUFFER_36, LABEL_STRING, "buffer_36");
        MDL::addLabel(BUFFER_37, LABEL_STRING, "buffer_37");
        MDL::addLabel(BUFFER_38, LABEL_STRING, "buffer_38");
        MDL::addLabel(BUFFER_39, LABEL_STRING, "buffer_39");
        MDL::addLabel(BUFFER_40, LABEL_STRING, "buffer_40");
        MDL::addLabel(BUFFER_41, LABEL_STRING, "buffer_41");
        MDL::addLabel(BUFFER_42, LABEL_STRING, "buffer_42");
        MDL::addLabel(BUFFER_43, LABEL_STRING, "buffer_43");
        MDL::addLabel(BUFFER_44, LABEL_STRING, "buffer_44");
        MDL::addLabel(BUFFER_45, LABEL_STRING, "buffer_45");
        MDL::addLabel(BUFFER_46, LABEL_STRING, "buffer_46");
        MDL::addLabel(BUFFER_47, LABEL_STRING, "buffer_47");
        MDL::addLabel(BUFFER_48, LABEL_STRING, "buffer_48");
        MDL::addLabel(BUFFER_49, LABEL_STRING, "buffer_49");
        MDL::addLabel(BUFFER_50, LABEL_STRING, "buffer_50");
        MDL::addLabel(BUFFER_51, LABEL_STRING, "buffer_51");
        MDL::addLabel(BUFFER_52, LABEL_STRING, "buffer_52");
        MDL::addLabel(BUFFER_53, LABEL_STRING, "buffer_53");
        MDL::addLabel(BUFFER_54, LABEL_STRING, "buffer_54");
        MDL::addLabel(BUFFER_55, LABEL_STRING, "buffer_55");
        MDL::addLabel(BUFFER_56, LABEL_STRING, "buffer_56");
        MDL::addLabel(BUFFER_57, LABEL_STRING, "buffer_57");
        MDL::addLabel(BUFFER_58, LABEL_STRING, "buffer_58");
        MDL::addLabel(BUFFER_59, LABEL_STRING, "buffer_59");
        MDL::addLabel(BUFFER_60, LABEL_STRING, "buffer_60");
        MDL::addLabel(BUFFER_61, LABEL_STRING, "buffer_61");
        MDL::addLabel(BUFFER_62, LABEL_STRING, "buffer_62");
        MDL::addLabel(BUFFER_63, LABEL_STRING, "buffer_63");
        MDL::addLabel(BUFFER_64, LABEL_STRING, "buffer_64");
        MDL::addLabel(BUFFER_65, LABEL_STRING, "buffer_65");
        MDL::addLabel(BUFFER_66, LABEL_STRING, "buffer_66");
        MDL::addLabel(BUFFER_67, LABEL_STRING, "buffer_67");
        MDL::addLabel(BUFFER_68, LABEL_STRING, "buffer_68");
        MDL::addLabel(BUFFER_69, LABEL_STRING, "buffer_69");
        MDL::addLabel(BUFFER_70, LABEL_STRING, "buffer_70");
        MDL::addLabel(BUFFER_71, LABEL_STRING, "buffer_71");
        MDL::addLabel(BUFFER_72, LABEL_STRING, "buffer_72");
        MDL::addLabel(BUFFER_73, LABEL_STRING, "buffer_73");
        MDL::addLabel(BUFFER_74, LABEL_STRING, "buffer_74");
        MDL::addLabel(BUFFER_75, LABEL_STRING, "buffer_75");
        MDL::addLabel(BUFFER_76, LABEL_STRING, "buffer_76");
        MDL::addLabel(BUFFER_77, LABEL_STRING, "buffer_77");
        MDL::addLabel(BUFFER_78, LABEL_STRING, "buffer_78");
        MDL::addLabel(BUFFER_79, LABEL_STRING, "buffer_79");
        MDL::addLabel(BUFFER_80, LABEL_STRING, "buffer_80");
        MDL::addLabel(BUFFER_81, LABEL_STRING, "buffer_81");
        MDL::addLabel(BUFFER_82, LABEL_STRING, "buffer_82");
        MDL::addLabel(BUFFER_83, LABEL_STRING, "buffer_83");
        MDL::addLabel(BUFFER_84, LABEL_STRING, "buffer_84");
        MDL::addLabel(BUFFER_85, LABEL_STRING, "buffer_85");
        MDL::addLabel(BUFFER_86, LABEL_STRING, "buffer_86");
        MDL::addLabel(BUFFER_87, LABEL_STRING, "buffer_87");
        MDL::addLabel(BUFFER_88, LABEL_STRING, "buffer_88");
        MDL::addLabel(BUFFER_89, LABEL_STRING, "buffer_89");
        MDL::addLabel(BUFFER_90, LABEL_STRING, "buffer_90");
        MDL::addLabel(BUFFER_91, LABEL_STRING, "buffer_91");
        MDL::addLabel(BUFFER_92, LABEL_STRING, "buffer_92");
        MDL::addLabel(BUFFER_93, LABEL_STRING, "buffer_93");
        MDL::addLabel(BUFFER_94, LABEL_STRING, "buffer_94");
        MDL::addLabel(BUFFER_95, LABEL_STRING, "buffer_95");
        MDL::addLabel(BUFFER_96, LABEL_STRING, "buffer_96");
        MDL::addLabel(BUFFER_97, LABEL_STRING, "buffer_97");
        MDL::addLabel(BUFFER_98, LABEL_STRING, "buffer_98");
        MDL::addLabel(BUFFER_99, LABEL_STRING, "buffer_99");

        // Add user-defined label aliases for allow use Xmipp-MetaData
        // with other labels
        MDL::addExtraAliases();
        MDL::bufferIndex = BUFFER_01;
    }

    ~MDLabelStaticInit()
    {
        //Free memory allocated for labels data
        for (int i = MDL_FIRST_LABEL; i < MDL_LAST_LABEL; ++i)
            delete MDL::data[i];
    }
    friend class MDL;
};
/** @} */

#endif
