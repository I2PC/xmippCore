/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *              Jan Horacek (xhorace4@fi.muni.cz)
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

#ifndef CORE_METADATALABEL_H
#define CORE_METADATALABEL_H

#include <map>
#include "xmipp_strings.h"
#include <vector>

class MDLabelData;
class MDObject;
class MDLabelStaticInit;

/** @addtogroup MetaData
 * @{
 */

/** Enumerate all posible labels to use in MetaData.
 */
enum MDLabel
{
    MDL_UNDEFINED = -1,
    MDL_FIRST_LABEL, ///< The label MDL_OBJID is special and should not be used
    MDL_OBJID = MDL_FIRST_LABEL, ///< object id (int), NOTE: This label is special and shouldn't be used
    MDL_GATHER_ID, /// Special label to be used when gathering MDs in MpiMetadataPrograms

    MDL_ANGLE_PSI, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI2, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI3, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI_DIFF, ///< difference between psi angles (double,degrees)
    MDL_ANGLE_ROT, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT2, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT3, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT_DIFF, ///< difference between rot angles (double,degrees)
    MDL_ANGLE_TILT, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT2, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT3, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT_DIFF, ///< difference between tilt angles (double,degrees)
    MDL_ANGLE_DIFF0, ///< difference between two angles (double,degrees)
    MDL_ANGLE_DIFF, ///< difference between two angles (double,degrees)
    MDL_ANGLE_DIFF2, ///< difference between two angles (double,degrees)
    MDL_ANGLE_Y,   ///< Angle between y-axis and tilt-axis (double, degrees) for untilted micrographs
    MDL_ANGLE_Y2,   ///< Angle between y-axis and tilt-axis (double, degrees) for tilted micrographs
    MDL_ANGLE_TEMPERATURE, ///< Angular temperature (double,degrees)
    MDL_APPLY_SHIFT,///<Apply shift when project the volume ,
    MDL_AVG, ///< average value (double)
    MDL_AVG_CHANGES_ORIENTATIONS, /// Average change in angular orientation (double degrees)
    MDL_AVG_CHANGES_OFFSETS, /// Average change in offset (double pixels)
    MDL_AVG_CHANGES_CLASSES, /// Average change in class assignment(double dimensionaless)
    MDL_AVGPMAX, ///< Average (per class) of the maximum value of normalized probability function) (double)
    MDL_BFACTOR, /// <Bfactor of a map, or even a local bfactor
    MDL_BGMEAN, ///< Mean background value for an image
    MDL_BLOCK_NUMBER, ///< Current block number (for incremental EM)

    MDL_CL2D_CHANGES, ///< Number of changes between iterations
    MDL_CL2D_SIMILARITY, ///< Average cross-correlation for the image (double)
    MDL_CLASS_COUNT, ///< Number of images assigned to the same class as this image
    MDL_CLASS_PERCENTAGE, ///< Percentage of images assigned to the same class as this image
    MDL_CLASSIFICATION_DATA, ///< Data vector for classification (vector double)
    MDL_CLASSIFICATION_DATA_SIZE, ///< Size of data vectors for classification (int)
    MDL_CLASSIFICATION_DPR_05, ///< Differential Phase Residual evaluated at FRC=0.5
    MDL_CLASSIFICATION_INTRACLASS_DISTANCE, ///< Average intraclass distance (double)
    MDL_CLASSIFICATION_FRC_05, ///< Digital frequency at which the FRC drops below 0.5 (double)
    MDL_COMMENT, ///< Serve to make annotations on the metadata row
    MDL_COORD_CONSENSUS_SCORE, ///< Store a score for the coords. consensus (it will change the behavoir of the viewer)
    MDL_COST, ///< Cost for the image (double)
    MDL_COST_PERCENTILE, ///< Cost percentile for the image (double)
    MDL_COUNT, ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]
    MDL_COUNT2, ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]
    MDL_CORR_DENOISED_PROJECTION, ///<Correlation between the denoised image and the projection proposed
    MDL_CORR_DENOISED_NOISY, ///<Correlation between the denoised image and the noisy version

    MDL_CLASS_INTERSECTION_SIZE_PVALUE, ///< P-value (1-percentile) of the class size compared to a random distribution (double in [0, 1])
    MDL_CLASS_INTERSECTION_RELATIVE_SIZE_PVALUE, ///< P-value (1-percentile) of the relative class size (size/max(original_classes)) compared to a random distribution (double in [0, 1])

    MDL_CRYSTAL_CELLX, ///< Cell location for crystals
    MDL_CRYSTAL_CELLY, ///< Cell location for crystals
    MDL_CRYSTAL_LATTICE_A,   /// < Lattice vector for projection (vector double)
    MDL_CRYSTAL_LATTICE_B,   /// < Lattice vector for projection (vector double)
    MDL_CRYSTAL_DISAPPEAR_THRE,   /// < Disappearing threshold (double)
    MDL_CRYSTAL_SHFILE,   /// < Shift file for crystal projection
    MDL_CRYSTAL_ORTHO_PRJ,   /// <Orthogonal projection or not (bool)
    MDL_CRYSTAL_PROJ,   /// < Have a crystal projection (bool)
    MDL_CRYSTAL_SHIFTX, ///< Shift for the image in the X axis (double) for crystals
    MDL_CRYSTAL_SHIFTY, ///< Shift for the image in the Y axis (double) for crystals
    MDL_CRYSTAL_SHIFTZ, ///< Shift for the image in the Z axis (double) for crystals
    MDL_CRYSTAL_NOISE_SHIFT , ///< noise if center of unit cell (vector double)

    MDL_CONTINUOUS_X, ///< X shift of continuous assignment
    MDL_CONTINUOUS_Y, ///< Y shift of continuous assignment
    MDL_CONTINUOUS_FLIP, ///< Flip of continuous assignment
    MDL_CONTINUOUS_GRAY_A, ///< a value of continuous assignment
    MDL_CONTINUOUS_GRAY_B, ///< b value of continuous assignment
    MDL_CONTINUOUS_SCALE_ANGLE, ///< scale angle of continuous assignment
    MDL_CONTINUOUS_SCALE_X, ///< scale x of continuous assignment
    MDL_CONTINUOUS_SCALE_Y, ///< scale y of continuous assignment

    MDL_CORRELATION_IDX, ///< correlation value between a particle and its assigned projection
    MDL_CORRELATION_MASK, ///< masked correlation value between a particle and its assigned projection inside the region with pixel values higher than the standard deviation
    MDL_CORRELATION_WEIGHT, ///< weighted correlation value between a particle and its assigned projection taking into the difference between both images

    MDL_CTF_DATA_PHASE_FLIPPED, // Is the Data Phase-Flippled?
    MDL_CTF_CORRECTED, // Is the CTF corrected?
    MDL_CTF_INPUTPARAMS, ///< Parameters file for the CTF Model (std::string)
    MDL_CTF_MODEL, ///< Name for the CTF Model (std::string)
    MDL_CTF_MODEL2, ///< Name for another CTF model (std::string)
    MDL_CTF_SAMPLING_RATE, ///< Sampling rate
    MDL_CTF_SAMPLING_RATE_Z, ///< Sampling rate in Z direction
    MDL_CTF_VOLTAGE, ///< Microscope voltage (kV)
    MDL_CTF_DEFOCUSA, ///< average defocus (Angtroms)
    MDL_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
    MDL_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
    MDL_CTF_DEFOCUS_CHANGE, ///< Defocus change with respect to previous defoucs (Angstroms)
    MDL_CTF_DEFOCUS_R2, ///< Defocus coefficient of determination
    MDL_CTF_DEFOCUS_COEFS, ///< Coefficients of the defocus adjustment plane
    MDL_CTF_DEFOCUS_RESIDUAL, ///< Difference between the observed defocus value and the estimated defocus value
    MDL_CTF_X0, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_Y0, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_XF, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_YF, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_DEFOCUS_PLANEUA, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEUB, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEUC, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEVA, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEVB, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEVC, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_ANGLE, ///< Defocus angle (degrees)
    MDL_CTF_CS, ///< Spherical aberration
    MDL_CTF_CA, ///< Chromatic aberration
    MDL_CTF_GROUP, ///< group images by defocus
    MDL_CTF_ENERGY_LOSS, ///< Energy loss
    MDL_CTF_ENVELOPE, //<Envelope function
    MDL_CTF_ENVELOPE_PLOT, //<Envelope function
    MDL_CTF_LENS_STABILITY, ///< Lens stability
    MDL_CTF_CONVERGENCE_CONE, ///< Convergence cone
    MDL_CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
    MDL_CTF_TRANSVERSAL_DISPLACEMENT, ///< Transversal displacemente
    MDL_CTF_Q0, ///< Inelastic absorption
    MDL_CTF_K, ///< CTF gain
    MDL_CTF_ENV_R0, ///< CTF Envelope polynomial parameter
    MDL_CTF_ENV_R1, ///< CTF Envelope polynomial parameter
    MDL_CTF_ENV_R2, ///< CTF Envelope polynomial parameter
    MDL_CTF_BG_GAUSSIAN_K, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_SIGMAU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_SIGMAV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_CU, ///<CTF_ CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_CV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_ANGLE, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_K, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_U, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_V, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_ANGLE, ///< CTF Background parameter
    MDL_CTF_BG_BASELINE, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_K, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_SIGMAU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_SIGMAV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_CU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_CV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_ANGLE, ///< CTF Background parameter
    MDL_CTF_BG_R1, ///< CTF Background polynomial parameter
    MDL_CTF_BG_R2, ///< CTF Background polynomial parameter
    MDL_CTF_BG_R3, ///< CTF Background polynomial parameter
    MDL_CTF_CRIT_NONASTIGMATICVALIDITY, ///< Maximum frequency (in Angstroms) at which non-astigmatic CTF correction is valid
    MDL_CTF_CRIT_PSDCORRELATION90, ///< PSD correlation at 90 degrees
    MDL_CTF_CRIT_FIRSTZERORATIO, ///< First zero ratio
    MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_RATIO, ///< Ratio sigma(firstMinimum)/sigma(firstZero)
    MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO, ///< Ratio sigma(firstMinimum-firstZero)/sigma(firstZero)
    MDL_CTF_CRIT_FIRSTZEROAVG, ///< First zero average (in Angstroms)
    MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, ///< First zero disagreement with second model (in Angstroms)
    MDL_CTF_CRIT_MAXFREQ, ///< Maximum frequency at which the envelope drops below 0.1 (in Angstroms)
    MDL_CTF_CRIT_DAMPING, ///< Minimum damping at border
    MDL_CTF_CRIT_PSDRADIALINTEGRAL, ///< Integral of the radial PSD
    MDL_CTF_CRIT_FITTINGSCORE, ///< Score of the fitting
    MDL_CTF_CRIT_FITTINGCORR13, ///< Correlation between the 1st and 3rd ring of the CTF
    MDL_CTF_CRIT_PSDVARIANCE, ///< PSD variance
    MDL_CTF_CRIT_PSDPCA1VARIANCE, ///< Variance in the first principal component of the PSDs
    MDL_CTF_CRIT_PSDPCARUNSTEST, ///< Runs test on the projection of the PSD on the first principal component
    MDL_CTF_CRIT_NORMALITY, ///< Normality test between histogram of micrography and gaussian distribution
    MDL_CTF_CRIT_ICENESS, ///< Iceness of the micrograph
    MDL_CTF_DOWNSAMPLE_PERFORMED, ///< Downsampling performed to estimate the CTF
    MDL_CTF_DIMENSIONS, // Size in pixels of the 3D PSF to be created (Xdim, Ydim, Zdim)
    MDL_CTF_LAMBDA, /// Wavelength (nm)
    MDL_CTF_XRAY_LENS_TYPE, ///Algorithm used to generate Xray PSF
    MDL_CTF_XRAY_OUTER_ZONE_WIDTH, /// Outermost zone width of the X-ray Fresnel lens (nm)
    MDL_CTF_XRAY_ZONES_NUMBER, // Number of zones of the X-ray Fresnel lens
    MDL_CTF_PHASE_SHIFT,    //Volta Phase Plate phase shift
    MDL_CTF_VPP_RADIUS,    //Phase Plate radius
    MDL_CUMULATIVE_SSNR, ///<Cumulative SSNR (double)
    MDL_DATATYPE, ///< if read from file original image datatype, this is an struct defined in image
    MDL_DEFGROUP, ///< Defocus group
    MDL_DIMRED, ///< Projection onto a reduced manifold (vector double)
    MDL_DIRECTION, ///< Direction in 3D


    MDL_DM3_IDTAG,
    MDL_DM3_NODEID,
    MDL_DM3_NUMBER_TYPE,
    MDL_DM3_PARENTID,
    MDL_DM3_TAGCLASS,
    MDL_DM3_TAGNAME,
    MDL_DM3_SIZE,
    MDL_DM3_VALUE,
    MDL_DOSE, ///< Dose e/A^2 (double)

    //End of labels

    MDL_ENABLED, ///< Is this image enabled? (int [-1 or 1])

    MDL_DATE,// < timestamp (string)
    MDL_TIME,// <  time in seconds (double)

    MDL_FLIP, ///< Flip the image? (bool)
    MDL_FOM, ///< Figure of Merit in 0-1 range (double)
    MDL_FRAME_ID, ///< Unique id of frame inside a Movie
	MDL_GRAPH_DISTANCE2MAX, ///< Distance to graph filtered max
	MDL_GRAPH_DISTANCE2MAX_PREVIOUS, ///< when previous assignment validation
	MDL_GRAPH_CC, ///< Correlation between assigned direction and graph filtered maximum
	MDL_GRAPH_CC_PREVIOUS, ///< when previous assignment validation
    MDL_HALF1, ///< Half map/image 1 file name (string)
    MDL_HALF2, ///< Half map/image file name (string)
    MDL_IDX, ///< Index within a list (size_t)
    MDL_IMAGE, ///< Name of an image (std::string)
    MDL_IMAGE_COVARIANCE, ///< Name of the covariance imagee associated to this image
    MDL_IMAGE_IDX, ///< Index of an image within a list (size_t)
    MDL_IMAGE_ORIGINAL, ///< Name of an image from which MDL_IMAGE is coming from
    MDL_IMAGE_REF, ///< Name of of the class image from which MDL_IMAGE is coming from
    MDL_IMAGE_RESIDUAL, ///< Name of a residual image associated to this image
    MDL_IMAGE_TILTED, ///< Name of the tilted images associated to MDL_IMAGE
    MDL_IMED, ///< imed value between a particle and its assigned projection
    MDL_IMGMD, ///< Name of Metadata file for all images (string)
    MDL_IMAGE1, ///< Image associated to this object (std::string)
    MDL_IMAGE2, ///< Image associated to this object (std::string)
    MDL_IMAGE3, ///< Image associated to this object (std::string)
    MDL_IMAGE4, ///< Image associated to this object (std::string)
    MDL_IMAGE5, ///< Image associated to this object (std::string)
    MDL_INTSCALE, ///< Intensity scale for an image
    MDL_ITEM_ID, ///< Unique identifier for items inside a list or set (std::size_t)
    MDL_ITER, ///< Current iteration number (int)
    MDL_KERDENSOM_FUNCTIONAL, ///< Functional value (double)
    MDL_KERDENSOM_REGULARIZATION, ///< Regularization value (double)
    MDL_KERDENSOM_SIGMA, ///< Sigma value (double)
    MDL_KEYWORDS, ///< Keywords associated to this line, should be a single string block (do not use spaces as separators)
    MDL_KMEANS2D_CENTROID, ///< Centroid of a cluster for the KMEANS2D classification
    MDL_KSTEST, ///<KS-test statistics
    MDL_LL, ///< contribution of an image to log-likelihood value
    MDL_LOCAL_ALIGNMENT_PATCHES, ///< Two values representing number of patches used for local alignment (X, Y)
    MDL_LOCAL_ALIGNMENT_COEFFS_X, ///< BSpline coefficient in X dim
    MDL_LOCAL_ALIGNMENT_COEFFS_Y, ///< BSpline coefficient in Y dim
    MDL_LOCAL_ALIGNMENT_CONF_2_5_PERC, ///< A shift amount at confidence level of 2.5%
    MDL_LOCAL_ALIGNMENT_CONF_97_5_PERC, ///< A shift amount at confidence level of 95.5%
    MDL_LOCAL_ALIGNMENT_CONTROL_POINTS, ///< Three values representing number of control points used for local alignment (X, Y, N)
    MDL_LOCAL_AVERAGE, ///< average value in the micrograph (double)
    MDL_MAGNIFICATION, /// Magnification of microscope
    MDL_MAPTOPOLOGY, ///< Map topology (KerDenSOM, ...)
    MDL_MASK, ///< Name of a mask associated to image
    MDL_MAXCC, ///< Maximum cross-correlation for the image (double)
    MDL_MAXCC_PERCENTILE, ///< Percentile of the maximum cross-correlation for the image (double)
    MDL_MAX, ///< Maximum value (double)
	MDL_MAXCC_PREVIOUS, ///< Correlation from previous alignment
    MDL_MICROGRAPH, ///< Name of a micrograph (std::string)
    MDL_MICROGRAPH_ID, ///< Micrograph unique id for reference (MDL_ITEM_ID should be used for Micrographs list)
    MDL_MICROGRAPH_MOVIE, ///< Name of a movie (std::string)
    MDL_MICROGRAPH_MOVIE_ID, ///< Unique identifier of a movie.
    MDL_MICROGRAPH_PARTICLES, ///< Name of a position file (std::string)
    MDL_MICROGRAPH_ORIGINAL, ///< Name of the original micrograph, MDL_MICROGRAPH is probably a downsampled version of this one (std::string)
    MDL_MICROGRAPH_TILTED, ///< Name of the corresponding tilted micrograph (std::string)
    MDL_MICROGRAPH_TILTED_ORIGINAL, ///< Name of the corresponding original, tilted micrograph (std::string)
    MDL_MIN, ///< Minimum value (double)
    MDL_MIRRORFRAC, ///< Mirror fraction for a Maximum Likelihood model
    MDL_MISSINGREGION_NR, ///< Number of missing region in subtomogram
    MDL_MISSINGREGION_TYPE, ///< Type of missing region in subtomogram
    MDL_MISSINGREGION_THY0, ///< Initial tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THYF, ///< Final tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THX0, ///< Initial tilt angle in X for missing region in subtomogram
    MDL_MISSINGREGION_THXF, ///< Final tilt angle in X for missing region in subtomogram

    MDL_MLF_CTF,    ///< MLF CTF estimated value (double)
    MDL_MLF_WIENER, ///< MLF Wiener filter (double)
    MDL_MLF_SIGNAL, ///< MLF signal (double)
    MDL_MLF_NOISE,  ///< MLF Wiener filter (double)

    MDL_MODELFRAC, ///< Model fraction (alpha_k) for a Maximum Likelihood model
    MDL_NEIGHBORS, ///< Vector of indexes to points some "neighbors"
    MDL_NEIGHBOR, ///< Particular neighbor (pointed myNEIGHBORS)
    MDL_NEIGHBORHOOD_RADIUS, ///< Radius of the neighborhood (radians)
    MDL_NMA, ///< Normal mode displacements (vector double)
    MDL_NMA_COLLECTIVITY, ///< NMA Collectivity of a given mode
    MDL_NMA_ATOMSHIFT, ///< NMA Atom shift in Angstroms
    MDL_NMA_ENERGY, ///< NMA energy contained in the NMA displacement vector
    MDL_NMA_MINRANGE, ///< Minimum value observed for a given NMA mode
    MDL_NMA_MAXRANGE, ///< Maximum value observed for a given NMA mode
    MDL_NMA_MODEFILE, ///< File with an NMA mode
    MDL_NMA_SCORE, ///< NMA Score of a given mode
    MDL_NMA_EIGENVAL, ///< NMA Eigenvalue of a given mode
    MDL_NOISE_ANGLES, ///< Noise description for projected angles
    MDL_NOISE_PARTICLE_COORD, ///< Noise description for particle's center coordenates (when projecting)
    MDL_NOISE_COORD,  //Use instead of MDL_NOISE_PARTICLE_COORD in future
    MDL_NOISE_PIXEL_LEVEL, ///< Noise description for pixels' gray level (when projecting)
    MDL_ORDER, /// auxiliary label to be used as an index (long)
    MDL_ORIGIN_X, ///< Origin for the image in the X axis (double)
    MDL_ORIGIN_Y, ///< Origin for the image in the Y axis (double)
    MDL_ORIGIN_Z, ///< Origin for the image in the Z axis (double)

    MDL_OPTICALFLOW_MEANX, ///<Mean of the movement in x direction of the motion map> (double)
    MDL_OPTICALFLOW_MEANY, ///<Mean of the movement in y direction of the motion map> (double)
    MDL_OPTICALFLOW_STDX, ///<Standatd deviation of the movement in x direction of the motion map> (double)
    MDL_OPTICALFLOW_STDY, ///<Standard deviation of the movement in y direction of the motion map> (double)

    MDL_PARTICLE_ID, ///< Particle unique identifier for reference. (The MDL_ITEM_ID should be used when particle list)
    MDL_PHANTOM_BGDENSITY, ///< Phantom background density (double)
    MDL_PHANTOM_FEATURE_CENTER, ///< Center of the feature (vector double)
    MDL_PHANTOM_FEATURE_DENSITY, ///< The density of the feature (double)
    MDL_PHANTOM_FEATURE_OPERATION, ///< Operation in case of overlapping features (+,-)
    MDL_PHANTOM_FEATURE_SPECIFIC, ///< Specific parameters for a feature (vector double)
    MDL_PHANTOM_FEATURE_TYPE, ///< Type of the feature (Sphere, Blob, ...) (std::string)
    MDL_PHANTOM_SCALE, ///< Number which will multiply all features (double)

    MDL_MACRO_CMD, //ImageJ macro command on picker
    MDL_MACRO_CMD_ARGS, //ImageJ macro args on picker
    MDL_COLOR, ///< Color for particle picking
    MDL_PICKING_TEMPLATES, ///< Number of templates
    MDL_PICKING_STATE, ///< State for particle picking
    MDL_PICKING_MICROGRAPH_STATE, ///< Micrograph state for particle picking
    MDL_PICKING_AUTOPICKPERCENT,
    MDL_PICKING_PARTICLE_SIZE, ///< Particle size for particle picking
    MDL_PICKING_AUTOPARTICLES_SIZE, ///< Number of automatic particles picked
    MDL_PICKING_MANUALPARTICLES_SIZE, ///< Number of manual particles picked
    MDL_PMAX, ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
    MDL_POINTSASYMETRICUNIT, /// < Number of non-redundant projections directions (size_t)

    MDL_PRJ_DIMENSIONS, // X,Y dimensions for the generated projections
    MDL_PRJ_ANGFILE,  ///< File for generated angles
    MDL_PRJ_PSI_NOISE,  /// < Psi angle dev and mean noise (vector double)
    MDL_PRJ_PSI_RANDSTR, /// < Type of randomness for Psi (std::string)
    MDL_PRJ_PSI_RANGE,  /// < Psi angle range (vector double)
    MDL_PRJ_ROT_NOISE ,  /// < Rotational angle dev and mean noise (vector double)
    MDL_PRJ_ROT_RANDSTR,  /// < Type of randomness for Rotational (std::string)
    MDL_PRJ_ROT_RANGE,
    MDL_PRJ_TILT_NOISE,  /// < Tilt angle dev and mean noise (vector double)
    MDL_PRJ_TILT_RANDSTR,  /// < Type of randomness for Tilt (std::string)
    MDL_PRJ_TILT_RANGE, // Vector with the initial and final tilt angle values, and step size
    MDL_PRJ_VOL,        // Volume file name to generate projections from

    MDL_PROGRAM,// <  program name
    MDL_USER,// <  user name

    MDL_DIMENSIONS_3D,  // X,Y,Z dimensions
    MDL_DIMENSIONS_2D,  // X,Y dimensions
    MDL_PSD, ///< A Power Spectrum Density file name (std::string)
    MDL_PSD_ENHANCED, ///< A enhanced Power Spectrum Density file name (std::string)
    MDL_RANDOMSEED, ///< Seed for random number generator
    MDL_REF3D, ///< 3D Class to which the image belongs (int)
    MDL_REF, ///< Class to which the image belongs (int)
    MDL_REF2, ///< Store a second class (int)
    MDL_REFMD, ///< Name of Metadata file for all references(string)
	MDL_ASSIGNED_DIR_REF_CC, ///< correlation of references assigned by two methods

    MDL_RESIDUE,        //<residue of an atomic model (int)
    MDL_RESOLUTION_ANISOTROPY, ///<Resolution anisotropy used to store the significance of the Bingham Test (double)
    MDL_RESOLUTION_DPR, ///<differential phase residual (double)
    MDL_RESOLUTION_ERRORL2, ///<Error in l2 (double)
    MDL_RESOLUTION_FRC, ///<Fourier shell correlation (double)
    MDL_RESOLUTION_FRCRANDOMNOISE, ///<Fourier shell correlation noise (double)
    MDL_RESOLUTION_FREQ, ///<Frequency in 1/A (double)
    MDL_RESOLUTION_FREQ2, ///< Frequency in 1/A squared (double)
    MDL_RESOLUTION_FREQREAL, ///< Frequency in A (double)
    MDL_RESOLUTION_FSO, ///<Fourier shell occupancy (double)
    MDL_RESOLUTION_LOCAL_RESIDUE, ///< Frequency in A of the local resolution of a residue (double)
    MDL_RESOLUTION_LOG_STRUCTURE_FACTOR, ///<Logarithm of the structure factor
    MDL_RESOLUTION_SSNR, ///<Fourier shell correlation (double)
    MDL_RESOLUTION_STRUCTURE_FACTOR, ///<Structure factor
    MDL_RESOLUTION_RFACTOR, /// rfactor

    MDL_SAMPLINGRATE, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_ORIGINAL, ///< original sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_X, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_Y, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_Z, ///< sampling rate in A/pixel (double)

    MDL_SCALE, ///< scaling factor for an image or volume (double)
    MDL_SCORE_BY_PCA_RESIDUAL_PROJ,
    MDL_SCORE_BY_PCA_RESIDUAL_EXP,
    MDL_SCORE_BY_PCA_RESIDUAL,
    MDL_SCORE_BY_ALIGNABILITY, ///< score by alignability (double)
    MDL_SCORE_BY_ALIGNABILITY_PRECISION, ///< score by alignability (double)
    MDL_SCORE_BY_ALIGNABILITY_ACCURACY, ///< score by alignability (double)
    MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP, ///< score by alignability experimental particles (double)
    MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF, ///< score by alignability references (double)
    MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP, ///< score by alignability experimental particles (double)
    MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF, ///< score by alignability references (double)
    MDL_SCORE_BY_ALIGNABILITY_NOISE, ///< score by alignability noise (double)
    MDL_SCORE_BY_EMPTINESS, ///< Small values represent worse particles. Much larger than 1 for good particles
    MDL_SCORE_BY_ENTROPY,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_GRANULO,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_HISTDIST,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_LBP,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_MIRROR, ///< score by mirror (double)
    MDL_SCORE_BY_RAMP,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_SCREENING, ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_VARIANCE,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_VAR, /// < Particle variance (double)
    MDL_SCORE_BY_GINI, /// < Micrographs Gini Coeff. (double)
    MDL_SCORE_BY_ZERNIKE,  ///< Feature vectors used to classify particles (vector double)
    MDL_SCORE_BY_ZSCORE,
    MDL_SELFILE, ///< Name of an image (std::string)
    MDL_SERIE, ///< A collection of micrographs, e.g. a tilt serie (std::string)
    MDL_SHIFT_X, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X2, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X3, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X_DIFF, ///< difference in Shift along X axis (double)
    MDL_SHIFT_Y, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y2, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y3, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y_DIFF, ///< difference in Shift along  Y axis (double)
    MDL_SHIFT_Z, ///< Shift for the image in the Z axis (double)
    MDL_SHIFT_Z2, ///< Shift for the image in the Z axis (double)
    MDL_SHIFT_Z3, ///< Shift for the image in the Z axis (double)
    MDL_SHIFT_DIFF0, ///< shift difference (double)
    MDL_SHIFT_DIFF, ///< shift difference (double)
    MDL_SHIFT_DIFF2, ///< shift difference (double)
    MDL_SIGMANOISE, ///< Standard deviation of the noise in ML model
    MDL_SIGMAOFFSET, ///< Standard deviation of the offsets in ML model
    MDL_SIGNALCHANGE, ///< Signal change for an image
    MDL_SPH_COEFFICIENTS, ///< Deformation coefficients
    MDL_SPH_DEFORMATION, ///< Deformation in voxels
    MDL_SPH_TSNE_COEFF1D, ///tsne coefficicient in 1D
    MDL_SPH_TSNE_COEFF2D, ///tsne coefficients in 2D
    MDL_STDDEV, ///<stdandard deviation value (double)
    MDL_STAR_COMMENT, ///< A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
	MDL_SUBTOMOID, ///<Subtomogram id (size_t)
    MDL_SUBTRACTION_R2, ///< R2 coefficient of subtracted particle 
    MDL_SUBTRACTION_B, ///< B (background) coefficient of adjusted model for subtract particle 
    MDL_SUBTRACTION_BETA0, ///< Beta 0 coefficient of adjusted model for subtract particle 
    MDL_SUBTRACTION_BETA1, ///< Beta 1 coefficient of adjusted model for subtract particle 
    MDL_SUM, ///< Sum of elements of a given type (double) [this is a genereic type do not use to transfer information to another program]
    MDL_SUMWEIGHT, ///< Sum of all weights in ML model
    MDL_SYMNO, ///< Symmetry number for a projection (used in ART)
	MDL_TILTPARTICLEID, ///<Tilt particle id (size_t)
    MDL_TOMOGRAM_VOLUME, ///< Name for the reconstructed tomogram volume (std::string)
    MDL_TOMOGRAMMD, ///< Name for a Metadata file (std::string)
    MDL_TRANSFORM_MATRIX, ///< transformation matrix in numpy string format or space separated (std::string)
    MDL_TSID, ///<Tilt series id (std::string)

    MDL_TEST_SIZE,// < number of test assigned to a program
    MDL_TILT_ANALYSIS_MEAN, // <Mean correlation between the PSD segments of a micrograph
    MDL_TILT_ANALYSIS_STD, // <STD correlation between the PSD segments of a micrograph
    MDL_TILT_ANALYSIS_MIN, // <Min correlation between the PSD segments of a micrograph
    MDL_TILT_ANALYSIS_MAX, // <Max correlation between the PSD segments of a micrograph
    MDL_TILT_ANALYSIS_PSDs, // <Image composed of the PSDs of each segment of a micrograph

    MDL_VOLUME_SCORE_SUM, /// < Score corresponding to the sum of cc with cc>threshold
    MDL_VOLUME_SCORE_SUM_TH, ///< Score corresponding to the sum of cc-threshold with cc>threshold
    MDL_VOLUME_SCORE_MEAN, ///< Score corresponding to the mean of cc with cc>threshold
    MDL_VOLUME_SCORE_MIN, ///< Score corresponding to the min of cc with cc>threshold
    MDL_VOLUME_SCORE1,/// < Score 1 for volumes
    MDL_VOLUME_SCORE2,/// < Score 2 for volumes
    MDL_VOLUME_SCORE3,/// < Score 3 for volumes
    MDL_VOLUME_SCORE4,/// < Score 4 for volumes
    MDL_WEIGHT, ///< Weight assigned to the image (double)
    MDL_WEIGHT_P, ///< Weight assigned to the image accordint to its clusterability with a significance with respect noise (double)
    MDL_WEIGHT_CONTINUOUS2, ///< Weight due to angular continuous assignment
    MDL_WEIGHT_JUMPER0, ///< Weight due to angular jumping
    MDL_WEIGHT_JUMPER, ///< Weight due to angular jumping
    MDL_WEIGHT_JUMPER2, ///< Weight due to angular jumping
    MDL_WEIGHT_REALCORR, ///< Weight due correlation between two images/subtomos
    MDL_WEIGHT_PHASECORR, ///< Weight due phase correlation between two images/subtomos
    MDL_WEIGHT_SIGNIFICANT, ///< Weight due to Angular significance
    MDL_WEIGHT_SSNR, ///< Weight due to SSNR
    MDL_WEIGHT_PRECISION_ALIGNABILITY, ///< Weight due to Alignability Precision
    MDL_WEIGHT_ALIGNABILITY, ///< Weight due to Alignability Precision and Accuracy
    MDL_WEIGHT_ACCURACY_ALIGNABILITY, ///< Weight due to Alignability Accuracy
    MDL_WEIGHT_PRECISION_MIRROR, ///< Weight due to Mirror Precision
    MDL_WROBUST, ///< Weight of t-student distribution in robust Maximum likelihood
    MDL_X, ///< X component (double)
    MDL_XCOOR, ///< X component (int)
    MDL_XCOOR_TILT, ///< X component in tilted micrograph (int)
    MDL_XSIZE, ///< X size (int)
    MDL_Y, ///< Y component (double)
    MDL_YCOOR, ///< Y component (int)
    MDL_YCOOR_TILT, ///< Y component in tilted micrograph (int)
    MDL_YSIZE, ///< Y size (int)
    MDL_Z, ///< Z component (double)
    MDL_ZCOOR, ///< Z component (int)
    MDL_ZSCORE, ///< Global Z Score (double)
    MDL_ZSCORE_DEEPLEARNING1, ///< Z Score (double)
    MDL_GOOD_REGION_SCORE, ///< Z Score (double)
    MDL_ZSCORE_HISTOGRAM, ///< Z Score (double)
    MDL_ZSCORE_RESMEAN, ///< Z Score of the mean of the residuals (double)
    MDL_ZSCORE_RESVAR, ///< Z Score of the stddev of the residuals (double)
    MDL_ZSCORE_RESCOV, ///< Z Score of the covariance matrix of the residuals (double)
    MDL_ZSCORE_SHAPE1, ///< Z Score (double)
    MDL_ZSCORE_SHAPE2, ///< Z Score (double)
    MDL_ZSCORE_SNR1, ///< Z Score (double)
    MDL_ZSCORE_SNR2, ///< Z Score (double)
    MDL_ZSIZE, ///< Z size (int)

    /** RELION labels */
    RLN_AREA_ID, ///< ID for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
    RLN_AREA_NAME, ///< Name for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
    RLN_COMMENT, // The RLN_COMMENT is handled specially as well

    RLN_CTF_BFACTOR, ///< B-factor
    RLN_CTF_SCALEFACTOR, ///< linear scale-factor
    RLN_CTF_SAMPLING_RATE, ///< Sampling rate
    RLN_CTF_VOLTAGE, ///< Microscope voltage (kV)
    RLN_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
    RLN_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
    RLN_CTF_DEFOCUS_ANGLE, ///< Defocus angle (degrees)
    RLN_CTF_CS, ///< Spherical aberration
    RLN_CTF_CA, ///< Chromatic aberration
    RLN_CTF_DETECTOR_PIXEL_SIZE, ///< Pixel size for detector as used in CTF-determination
    RLN_CTF_ENERGY_LOSS, ///< Energy loss
    RLN_CTF_FOM, ///< ctffind3 FOM (CC) for quality of CTF-fit
    RLN_CTF_IMAGE, ///< name of an image describing the CTF model
    RLN_CTF_LENS_STABILITY, ///< Lens stability
    RLN_CTF_MAGNIFICATION, ///< Magnification used for CTF-determination
    RLN_CTF_MAXRES, ///< Maximum resolution with Thon rings
    RLN_CTF_CONVERGENCE_CONE, ///< Convergence cone
    RLN_CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
    RLN_CTF_TRANSVERSAL_DISPLACEMENT, ///< Transversal displacemente
    RLN_CTF_Q0, ///< Amplitude contrast
    RLN_CTF_K, ///< CTF gain
    RLN_CTF_VALUE, ///< CTF value
    RLN_CTF_VALIDATIONSCORE, ///< Gctf-based validation score for CTF fit
    RLN_CTF_PHASESHIFT, ///< Phase-shift from a phase-plate (in degrees)

    RLN_IMAGE_NAME,
    RLN_IMAGE_RECONSTRUCT_NAME,
    RLN_IMAGE_ID,
    RLN_IMAGE_ENABLED,
    RLN_IMAGE_DATATYPE,
    RLN_IMAGE_DIMENSIONALITY,
    RLN_IMAGE_BEAMTILT_X,
    RLN_IMAGE_BEAMTILT_Y,
    RLN_IMAGE_BEAMTILT_GROUP,
    RLN_IMAGE_COORD_X,
    RLN_IMAGE_COORD_Y,
    RLN_IMAGE_COORD_Z,
    RLN_IMAGE_FRAME_NR,
    RLN_IMAGE_MAGNIFICATION_CORRECTION,
    RLN_IMAGE_NORM_CORRECTION,
    RLN_IMAGE_ORI_NAME,
    RLN_IMAGE_SAMPLINGRATE,
    RLN_IMAGE_SAMPLINGRATE_X,
    RLN_IMAGE_SAMPLINGRATE_Y,
    RLN_IMAGE_SAMPLINGRATE_Z,
    RLN_IMAGE_SIZE,
    RLN_IMAGE_SIZEX,
    RLN_IMAGE_SIZEY,
    RLN_IMAGE_SIZEZ,
    RLN_IMAGE_STATS_MIN,
    RLN_IMAGE_STATS_MAX,
    RLN_IMAGE_STATS_AVG,
    RLN_IMAGE_STATS_STDDEV,
    RLN_IMAGE_STATS_SKEW,
    RLN_IMAGE_STATS_KURT,
    RLN_IMAGE_WEIGHT,

    RLN_MASK_NAME,

    RLN_MATRIX_1_1,
    RLN_MATRIX_1_2,
    RLN_MATRIX_1_3,
    RLN_MATRIX_2_1,
    RLN_MATRIX_2_2,
    RLN_MATRIX_2_3,
    RLN_MATRIX_3_1,
    RLN_MATRIX_3_2,
    RLN_MATRIX_3_3,

    RLN_MICROGRAPH_ID,
    RLN_MICROGRAPH_MOVIE_NAME,
    RLN_MICROGRAPH_NAME,
    RLN_MICROGRAPH_NAME_WODOSE,
    RLN_MICROGRAPH_TILT_ANGLE,
    RLN_MICROGRAPH_TILT_AXIS_DIRECTION,
    RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE,

    RLN_MLMODEL_ACCURACY_ROT,
    RLN_MLMODEL_ACCURACY_TRANS,
    RLN_MLMODEL_AVE_PMAX,
    RLN_MLMODEL_CURRENT_RESOLUTION,
    RLN_MLMODEL_CURRENT_SIZE,
    RLN_MLMODEL_DATA_VS_PRIOR_REF,
    RLN_MLMODEL_DIMENSIONALITY,
    RLN_MLMODEL_DIMENSIONALITY_DATA,
    RLN_MLMODEL_DIFF2_HALVES_REF,
    RLN_MLMODEL_ESTIM_RESOL_REF,
    RLN_MLMODEL_FOURIER_COVERAGE_REF,
    RLN_MLMODEL_FOURIER_COVERAGE_TOTAL_REF,
    RLN_MLMODEL_FSC_HALVES_REF,
    RLN_MLMODEL_GROUP_NAME,
    RLN_MLMODEL_GROUP_NO,
    RLN_MLMODEL_GROUP_NR_PARTICLES,
    RLN_MLMODEL_GROUP_SCALE_CORRECTION,
    RLN_MLMODEL_HELICAL_NR_ASU,
    RLN_MLMODEL_HELICAL_TWIST,
    RLN_MLMODEL_HELICAL_TWIST_MIN,
    RLN_MLMODEL_HELICAL_TWIST_MAX,
    RLN_MLMODEL_HELICAL_TWIST_INITIAL_STEP,
    RLN_MLMODEL_HELICAL_RISE,
    RLN_MLMODEL_HELICAL_RISE_MIN,
    RLN_MLMODEL_HELICAL_RISE_MAX,
    RLN_MLMODEL_HELICAL_RISE_INITIAL_STEP,
    RLN_MLMODEL_IS_HELIX,
    RLN_MLMODEL_INTERPOLATOR,
    RLN_MLMODEL_LL,
    RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,
    RLN_MLMODEL_NORM_CORRECTION_AVG,
    RLN_MLMODEL_NR_BODIES,
    RLN_MLMODEL_NR_CLASSES,
    RLN_MLMODEL_NR_GROUPS,
    RLN_MLMODEL_ORIGINAL_SIZE,
    RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION,
    RLN_MLMODEL_PADDING_FACTOR,
    RLN_MLMODEL_PDF_CLASS,
    RLN_MLMODEL_PRIOR_OFFX_CLASS,
    RLN_MLMODEL_PRIOR_OFFY_CLASS,
    RLN_MLMODEL_PDF_ORIENT,
    RLN_MLMODEL_PIXEL_SIZE,
    RLN_MLMODEL_POWER_REF,
    RLN_MLMODEL_PRIOR_MODE,
    RLN_MLMODEL_SGD_GRADIENT_IMAGE,
    RLN_MLMODEL_SIGMA_OFFSET,
    RLN_MLMODEL_SIGMA_ROT,
    RLN_MLMODEL_SIGMA_TILT,
    RLN_MLMODEL_SIGMA_PSI,
    RLN_MLMODEL_REF_IMAGE,
    RLN_MLMODEL_SIGMA2_NOISE,
    RLN_MLMODEL_SIGMA2_REF,
    RLN_MLMODEL_SSNR_REF,
    RLN_MLMODEL_TAU2_FUDGE_FACTOR,
    RLN_MLMODEL_TAU2_REF,

    RLN_OPTIMISER_ACCURACY_ROT,
    RLN_OPTIMISER_ACCURACY_TRANS,
    RLN_OPTIMISER_ADAPTIVE_FRACTION,
    RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING,
    RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER,
    RLN_OPTIMISER_AVAILABLE_MEMORY,
    RLN_OPTIMISER_BEST_RESOL_THUS_FAR,
    RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
    RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS,
    RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES,
    RLN_OPTIMISER_COARSE_SIZE,
    RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED,
    RLN_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED,
    RLN_OPTIMISER_DATA_STARFILE,
    RLN_OPTIMISER_DO_AUTO_REFINE,
    RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES,
    RLN_OPTIMISER_DO_CORRECT_CTF,
    RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION,
    RLN_OPTIMISER_DO_CORRECT_NORM,
    RLN_OPTIMISER_DO_CORRECT_SCALE,
    RLN_OPTIMISER_DO_HELICAL_REFINE,
    RLN_OPTIMISER_DO_REALIGN_MOVIES,
    RLN_OPTIMISER_DO_MAP,
    RLN_OPTIMISER_DO_SGD,
    RLN_OPTIMISER_DO_SOLVENT_FLATTEN,
    RLN_OPTIMISER_DO_SKIP_ALIGN,
    RLN_OPTIMISER_DO_SKIP_ROTATE,
    RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES,
    RLN_OPTIMISER_DO_ZERO_MASK,
    RLN_OPTIMISER_FIX_SIGMA_NOISE,
    RLN_OPTIMISER_FIX_SIGMA_OFFSET,
    RLN_OPTIMISER_FIX_TAU,
    RLN_OPTIMISER_HAS_CONVERGED,
    RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT,
    RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO,
    RLN_OPTIMISER_HELICAL_TWIST_INITIAL,
    RLN_OPTIMISER_HELICAL_RISE_INITIAL,
    RLN_OPTIMISER_HELICAL_Z_PERCENTAGE,
    RLN_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER,
    RLN_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER,
    RLN_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT,
    RLN_OPTIMISER_HELICAL_SIGMA_DISTANCE,
    RLN_OPTIMISER_HIGHRES_LIMIT_SGD,
    RLN_OPTIMISER_IGNORE_HELICAL_SYMMETRY,
    RLN_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED,
    RLN_OPTIMISER_HIGHRES_LIMIT_EXP,
    RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK,
    RLN_OPTIMISER_INCR_SIZE,
    RLN_OPTIMISER_ITERATION_NO,
    RLN_OPTIMISER_LOCAL_SYMMETRY_FILENAME,
    RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES,
    RLN_OPTIMISER_MAGNIFICATION_RANGE,
    RLN_OPTIMISER_MAGNIFICATION_STEP,
    RLN_OPTIMISER_MAX_COARSE_SIZE,
    RLN_OPTIMISER_MAX_NR_POOL,
    RLN_OPTIMISER_MODEL_STARFILE,
    RLN_OPTIMISER_MODEL_STARFILE2,
    RLN_OPTIMISER_NR_ITERATIONS,
    RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN,
    RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES,
    RLN_OPTIMISER_OUTPUT_ROOTNAME,
    RLN_OPTIMISER_PARTICLE_DIAMETER,
    RLN_OPTIMISER_RADIUS_MASK_3D_MAP,
    RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES,
    RLN_OPTIMISER_RANDOM_SEED,
    RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED,
    RLN_OPTIMISER_SAMPLING_STARFILE,
    RLN_OPTIMISER_SGD_MU,
    RLN_OPTIMISER_SGD_SIGMA2FUDGE_INI,
    RLN_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE,
    RLN_OPTIMISER_SGD_SUBSET_START,
    RLN_OPTIMISER_SGD_SUBSET_SIZE,
    RLN_OPTIMISER_SGD_WRITE_EVERY_SUBSET,
    RLN_OPTIMISER_SGD_MAX_SUBSETS,
    RLN_OPTIMISER_SGD_STEPSIZE,
    RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES,
    RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS,
    RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS,
    RLN_OPTIMISER_SOLVENT_MASK_NAME,
    RLN_OPTIMISER_SOLVENT_MASK2_NAME,
    RLN_OPTIMISER_TAU_SPECTRUM_NAME,
    RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING,
    RLN_OPTIMISER_WIDTH_MASK_EDGE,

    RLN_ORIENT_FLIP,
    RLN_ORIENT_ID,
    RLN_ORIENT_ORIGIN_X,
    RLN_ORIENT_ORIGIN_X_PRIOR,
    RLN_ORIENT_ORIGIN_Y,
    RLN_ORIENT_ORIGIN_Y_PRIOR,
    RLN_ORIENT_ORIGIN_Z,
    RLN_ORIENT_ORIGIN_Z_PRIOR,
    RLN_ORIENT_ROT,
    RLN_ORIENT_ROT_PRIOR,
    RLN_ORIENT_TILT,
    RLN_ORIENT_TILT_PRIOR,
    RLN_ORIENT_PSI,
    RLN_ORIENT_PSI_PRIOR,
    RLN_ORIENT_PSI_PRIOR_FLIP_RATIO,

    RLN_PARTICLE_AUTOPICK_FOM,
    RLN_PARTICLE_CLASS,
    RLN_PARTICLE_DLL,
    RLN_PARTICLE_ID,
    RLN_PARTICLE_FOM,
    RLN_PARTICLE_KL_DIVERGENCE,
    RLN_PARTICLE_MOVIE_RUNNING_AVG,
    RLN_PARTICLE_RANDOM_SUBSET,
    RLN_PARTICLE_NAME,
    RLN_PARTICLE_ORI_NAME,
    RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES,
    RLN_PARTICLE_NR_FRAMES,
    RLN_PARTICLE_NR_FRAMES_AVG,

    RLN_PARTICLE_PMAX,

    RLN_PARTICLE_HELICAL_TUBE_ID,
    RLN_PARTICLE_HELICAL_TUBE_PITCH,
    RLN_PARTICLE_HELICAL_TRACK_LENGTH,

    RLN_PIPELINE_JOB_COUNTER,
    RLN_PIPELINE_NODE_NAME,
    RLN_PIPELINE_NODE_TYPE,
    RLN_PIPELINE_PROCESS_ALIAS,
    RLN_PIPELINE_PROCESS_NAME,
    RLN_PIPELINE_PROCESS_TYPE,
    RLN_PIPELINE_PROCESS_STATUS,
    RLN_PIPELINE_EDGE_FROM,
    RLN_PIPELINE_EDGE_TO,
    RLN_PIPELINE_EDGE_PROCESS,

    RLN_POSTPROCESS_AMPLCORR_MASKED,
    RLN_POSTPROCESS_AMPLCORR_UNMASKED,
    RLN_POSTPROCESS_BFACTOR,
    RLN_POSTPROCESS_DPR_MASKED,
    RLN_POSTPROCESS_DPR_UNMASKED,
    RLN_POSTPROCESS_FINAL_RESOLUTION,
    RLN_POSTPROCESS_FSC_GENERAL,
    RLN_POSTPROCESS_FSC_TRUE,
    RLN_POSTPROCESS_FSC_MASKED,
    RLN_POSTPROCESS_FSC_UNMASKED,
    RLN_POSTPROCESS_FSC_RANDOM_MASKED,
    RLN_POSTPROCESS_GUINIER_FIT_CORRELATION,
    RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT,
    RLN_POSTPROCESS_GUINIER_FIT_SLOPE,
    RLN_POSTPROCESS_GUINIER_VALUE_IN,
    RLN_POSTPROCESS_GUINIER_VALUE_INVMTF,
    RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED,
    RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED,
    RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT,
    RLN_POSTPROCESS_GUINIER_RESOL_SQUARED,
    RLN_POSTPROCESS_MTF_VALUE, ///< Detector MTF value

    RLN_SAMPLING_IS_3D,
    RLN_SAMPLING_IS_3D_TRANS,
    RLN_SAMPLING_HEALPIX_ORDER,
    RLN_SAMPLING_HELICAL_OFFSET_STEP,
    RLN_SAMPLING_LIMIT_TILT,
    RLN_SAMPLING_OFFSET_RANGE,
    RLN_SAMPLING_OFFSET_STEP,
    RLN_SAMPLING_PERTURB,
    RLN_SAMPLING_PERTURBATION_FACTOR,
    RLN_SAMPLING_PRIOR_MODE,
    RLN_SAMPLING_PSI_STEP,
    RLN_SAMPLING_SIGMA_ROT,
    RLN_SAMPLING_SIGMA_TILT,
    RLN_SAMPLING_SIGMA_PSI,
    RLN_SAMPLING_SYMMETRY,

    RLN_SELECTED,
    RLN_SELECT_PARTICLES_ZSCORE,
    RLN_SORTED_IDX,
    RLN_STARFILE_MOVIE_PARTICLES,
    RLN_PERFRAME_CUMULATIVE_WEIGHT,
    RLN_PERFRAME_RELATIVE_WEIGHT,

    RLN_RESOLUTION,
    RLN_RESOLUTION_ANGSTROM,
    RLN_RESOLUTION_INVPIXEL,
    RLN_SPECTRAL_IDX,

    /** Buffer labels to read generic labels */
    BUFFER_LABELS_START,
    BUFFER_01,
    BUFFER_02,
    BUFFER_03,
    BUFFER_04,
    BUFFER_05,
    BUFFER_06,
    BUFFER_07,
    BUFFER_08,
    BUFFER_09,
    BUFFER_10,
    BUFFER_11,
    BUFFER_12,
    BUFFER_13,
    BUFFER_14,
    BUFFER_15,
    BUFFER_16,
    BUFFER_17,
    BUFFER_18,
    BUFFER_19,
    BUFFER_20,
    BUFFER_21,
    BUFFER_22,
    BUFFER_23,
    BUFFER_24,
    BUFFER_25,
    BUFFER_26,
    BUFFER_27,
    BUFFER_28,
    BUFFER_29,
    BUFFER_30,
    BUFFER_31,
    BUFFER_32,
    BUFFER_33,
    BUFFER_34,
    BUFFER_35,
    BUFFER_36,
    BUFFER_37,
    BUFFER_38,
    BUFFER_39,
    BUFFER_40,
    BUFFER_41,
    BUFFER_42,
    BUFFER_43,
    BUFFER_44,
    BUFFER_45,
    BUFFER_46,
    BUFFER_47,
    BUFFER_48,
    BUFFER_49,
    BUFFER_50,
    BUFFER_51,
    BUFFER_52,
    BUFFER_53,
    BUFFER_54,
    BUFFER_55,
    BUFFER_56,
    BUFFER_57,
    BUFFER_58,
    BUFFER_59,
    BUFFER_60,
    BUFFER_61,
    BUFFER_62,
    BUFFER_63,
    BUFFER_64,
    BUFFER_65,
    BUFFER_66,
    BUFFER_67,
    BUFFER_68,
    BUFFER_69,
    BUFFER_70,
    BUFFER_71,
    BUFFER_72,
    BUFFER_73,
    BUFFER_74,
    BUFFER_75,
    BUFFER_76,
    BUFFER_77,
    BUFFER_78,
    BUFFER_79,
    BUFFER_80,
    BUFFER_81,
    BUFFER_82,
    BUFFER_83,
    BUFFER_84,
    BUFFER_85,
    BUFFER_86,
    BUFFER_87,
    BUFFER_88,
    BUFFER_89,
    BUFFER_90,
    BUFFER_91,
    BUFFER_92,
    BUFFER_93,
    BUFFER_94,
    BUFFER_95,
    BUFFER_96,
    BUFFER_97,
    BUFFER_98,
    BUFFER_99,

    MDL_LAST_LABEL  // **** NOTE ****: Do keep this label always at the end,it is here for looping purposes
};//close enum Label

typedef std::vector<MDLabel> MDLabelVector;

/** Macro for iterate over all labels */
#define FOR_ALL_LABELS() for (int _label = MDL_FIRST_LABEL; _label < MDL_LAST_LABEL; ++_label)

/** Possible types of the values of labels */
enum MDLabelType
{
    LABEL_NOTYPE = -1,
    LABEL_INT,
    LABEL_BOOL,
    LABEL_DOUBLE,
    LABEL_STRING,
    LABEL_VECTOR_DOUBLE,
    LABEL_SIZET,
    LABEL_VECTOR_SIZET,
    LABEL_VECTOR_FLOAT
};

/** Possible types of the values of labels */
enum MDLabelTag
{
    TAGLABEL_NOTAG = 0,
    TAGLABEL_TEXTFILE=0x1,
    TAGLABEL_METADATA=0x3,
    TAGLABEL_CTFPARAM=0x5,
    TAGLABEL_IMAGE=0x8,
    TAGLABEL_VOLUME=0x10,
    TAGLABEL_STACK=0x20,
    TAGLABEL_MICROGRAPH=0x48,
    TAGLABEL_PSD=0x88
};

/**Just an utility function */
bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label);

//Just an struct to store type and string alias
class MDLabelData
{
public:
    MDLabelType type;
    String str;
    int tags;
    //Default constructor
    MDLabelData()
    {
        type = LABEL_NOTYPE;
        tags=TAGLABEL_NOTAG;
    }

    MDLabelData(MDLabelType t, const String &s, int tags)
    {
        type = t;
        str = s;
        this->tags=tags;
    }
}
;//close class MDLabelData

/** Explicit instantiation */
#ifndef __APPLE__
template class std::vector<MDObject *>
;
#endif

#endif
