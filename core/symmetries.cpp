/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <stdio.h>

#include "symmetries.h"

// Read Symmetry file ======================================================
// crystal symmetry matices from http://cci.lbl.gov/asu_gallery/
int SymList::readSymmetryFile(FileName fn_sym, double accuracy)
{
    int i, j;
    FILE *fpoii;
    char line[80];
    char *auxstr;
    double ang_incr, rot_ang;
    int  fold;
    Matrix2D<double> L(4, 4), R(4, 4);
    Matrix1D<double> axis(3), shift(3);
    int pgGroup, pgOrder;
    std::vector<std::string> fileContent;

    //check if reserved word

    // Open file ---------------------------------------------------------
    if ((fpoii = fopen(fn_sym.c_str(), "r")) == NULL)
    {
        //check if reserved word and return group and order
        if (isSymmetryGroup(fn_sym, pgGroup, pgOrder))
        {
            fillSymmetryClass(fn_sym, pgGroup, pgOrder,fileContent);
        }
        else
            REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"SymList::read_sym_file:Can't open file: "
                         + " or do not recognize symmetry group" + fn_sym);
    }
    else
    {
        while (fgets(line, 79, fpoii) != NULL)
        {
            if (line[0] == ';' || line[0] == '#' || line[0] == '\0')
                continue;
            fileContent.push_back(line);
        }
        fclose(fpoii);
    }

    //reset space_group
    space_group = 0;
    // Count the number of symmetries ------------------------------------
    true_symNo = 0;
    // count number of axis and mirror planes. It will help to identify
    // the crystallographic symmetry

    int no_axis, no_mirror_planes, no_inversion_points;
    no_axis = no_mirror_planes = no_inversion_points = 0;

    for (size_t n=0; n<fileContent.size(); n++)
    {
        strcpy(line,fileContent[n].c_str());
        auxstr = firstToken(line);
        if (auxstr == NULL)
        {
            std::cout << line;
            std::cout << "Wrong line in symmetry file, the line is skipped\n";
            continue;
        }
        if (strcmp(auxstr, "rot_axis") == 0)
        {
            auxstr = nextToken();
            fold = textToInteger(auxstr);
            true_symNo += (fold - 1);
            no_axis++;
        }
        else if (strcmp(auxstr, "mirror_plane") == 0)
        {
            true_symNo++;
            no_mirror_planes++;
        }
        else if (strcmp(auxstr, "inversion") == 0)
        {
            true_symNo += 1;
            no_inversion_points = 1;
        }
        else if (strcmp(auxstr, "P4212") == 0)
            true_symNo += 7;
        else if (strcmp(auxstr, "P2_122") == 0)
            true_symNo += 3;
        else if (strcmp(auxstr, "P22_12") == 0)
            true_symNo += 3;
    }
    // Ask for memory
    __L.resize(4*true_symNo, 4);
    __R.resize(4*true_symNo, 4);
    __shift.resize(true_symNo, 3);
    __chain_length.resize(true_symNo);
    __chain_length.initConstant(1);

    // Read symmetry parameters
    i = 0;
    shift.initZeros();

    for (size_t n=0; n<fileContent.size(); n++)
    {
        strcpy(line,fileContent[n].c_str());
        auxstr = firstToken(line);
        // Rotational axis ---------------------------------------------------
        if (strcmp(auxstr, "rot_axis") == 0)
        {
            auxstr = nextToken();
            fold = textToInteger(auxstr);
            auxstr = nextToken();
            XX(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            YY(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            ZZ(axis) = textToFloat(auxstr);
            ang_incr = 360. / fold;
            L.initIdentity();
            for (j = 1, rot_ang = ang_incr; j < fold; j++, rot_ang += ang_incr)
            {
                rotation3DMatrix(rot_ang, axis, R);
                setShift(i, shift);
                setMatrices(i++, L, R.transpose());
            }
            __sym_elements++;
            // inversion ------------------------------------------------------
        }
        else if (strcmp(auxstr, "inversion") == 0)
        {
            L.initIdentity();
            L(2, 2) = -1;
            R.initIdentity();
            R(0, 0) = -1.;
            R(1, 1) = -1.;
            R(2, 2) = -1.;
            setShift(i, shift);
            setMatrices(i++, L, R);
            __sym_elements++;
            // mirror plane -------------------------------------------------------------
        }
        else if (strcmp(auxstr, "mirror_plane") == 0)
        {
            auxstr = nextToken();
            XX(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            YY(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            ZZ(axis) = textToFloat(auxstr);
            L.initIdentity();
            L(2, 2) = -1;
            Matrix2D<double> A;
            alignWithZ(axis,A);
            A = A.transpose();
            R = A * L * A.inv();
            setShift(i, shift);
            L.initIdentity();
            setMatrices(i++, L, R);
            __sym_elements++;
            // P4212 -------------------------------------------------------------
        }
        else if (strcmp(auxstr, "P4212") == 0)
        {
            space_group = sym_P42_12;
            accuracy = -1; // Do not compute subgroup
            L.initIdentity();

            // With 0 shift
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = R(1, 1) = -1;
            R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(2, 2) = -1;
            R(0, 1) = R(1, 0) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(2, 2) = R(0, 1) = R(1, 0) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            // With 1/2 shift
            VECTOR_R3(shift, 0.5, 0.5, 0);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 1) = -1;
            R(1, 0) = R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(1, 0) = -1;
            R(0, 1) = R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = R(2, 2) = -1;
            R(1, 1) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(1, 1) = R(2, 2) = -1;
            R(0, 0) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            __sym_elements++;
        }
        else if (strcmp(auxstr, "P2_122") == 0)
        {
            space_group = sym_P2_122;
            accuracy = -1; // Do not compute subgroup
            L.initIdentity();

            // With 0 shift
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = -1;
            R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            // With 1/2 shift
            VECTOR_R3(shift, 0.5, 0.0, 0.0);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = 1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = 1;
            R(1, 1) = -1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            __sym_elements++;
        }
        else if (strcmp(auxstr, "P22_12") == 0)
        {
            space_group = sym_P22_12;
            accuracy = -1; // Do not compute subgroup
            L.initIdentity();

            // With 0 shift
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = -1;
            R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            // With 1/2 shift
            VECTOR_R3(shift, 0.0, 0.5, 0.0);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = 1;
            R(1, 1) = -1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = 1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            __sym_elements++;
        }
    }

    if (accuracy > 0)
        computeSubgroup(accuracy);

    //possible crystallographic symmetry
    if (no_axis == 0 && no_mirror_planes == 0 && no_inversion_points == 0 &&
        true_symNo == 7 && space_group == sym_P42_12)
        space_group = sym_P42_12;
    else if (no_axis == 0 && no_mirror_planes == 0 && no_inversion_points == 0 &&
             true_symNo == 3 && space_group == sym_P2_122)
        space_group = sym_P2_122;
    else if (no_axis == 0 && no_mirror_planes == 0 && no_inversion_points == 0 &&
             true_symNo == 3 && space_group == sym_P22_12)
        space_group = sym_P22_12;
    // P4 and P6
    else if (no_axis == 1 && no_mirror_planes == 0 && no_inversion_points == 0 &&
             fabs(R(2, 2) - 1.) < XMIPP_EQUAL_ACCURACY &&
             fabs(R(0, 0) - R(1, 1)) < XMIPP_EQUAL_ACCURACY &&
             fabs(R(0, 1) + R(1, 0)) < XMIPP_EQUAL_ACCURACY)
    {
        switch (true_symNo)
        {
        case(5):
                        space_group = sym_P6;
            break;
        case(3):
                        space_group = sym_P4;
            break;
        default:
            space_group = sym_undefined;
            break;
        }//switch end
    }//end else if (no_axis==1 && no_mirror_planes== 0
    else if (no_axis == 0 && no_inversion_points == 0 && no_mirror_planes == 0)
        space_group = sym_P1;
    else
        space_group = sym_undefined;
    return pgGroup;
}

// Get matrix ==============================================================
void SymList::getMatrices(int i, Matrix2D<double> &L, Matrix2D<double> &R,
                          bool homogeneous)
const
{
    int k, kp, l;
    if (homogeneous)
    {
        L.initZeros(4, 4);
        R.initZeros(4, 4);
        for (k = 4 * i, kp=0; k < 4*i + 4; k++, kp++)
            for (l = 0; l < 4; l++)
            {
                dMij(L,kp, l) = dMij(__L,k, l);
                dMij(R,kp, l) = dMij(__R,k, l);
            }
    }
    else
    {
        L.initZeros(3, 3);
        R.initZeros(3, 3);
        for (k = 4 * i, kp=0; k < 4*i + 3; k++, kp++)
            for (l = 0; l < 3; l++)
            {
                dMij(L,kp, l) = dMij(__L,k, l);
                dMij(R,kp, l) = dMij(__R,k, l);
            }
    }
}

// Set matrix ==============================================================
void SymList::setMatrices(int i, const Matrix2D<double> &L,
                          const Matrix2D<double> &R)
{
    int k, l;
    for (k = 4 * i; k < 4*i + 4; k++)
        for (l = 0; l < 4; l++)
        {
            __L(k, l) = L(k - 4 * i, l);
            __R(k, l) = R(k - 4 * i, l);
        }
}

// Get/Set shift ===========================================================
void SymList::getShift(int i, Matrix1D<double> &shift) const
{
    shift.resize(3);
    XX(shift) = __shift(i, 0);
    YY(shift) = __shift(i, 1);
    ZZ(shift) = __shift(i, 2);
}

void SymList::setShift(int i, const Matrix1D<double> &shift)
{
    if (shift.size() != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "SymList::add_shift: Shift vector is not 3x1");
    __shift(i, 0) = XX(shift);
    __shift(i, 1) = YY(shift);
    __shift(i, 2) = ZZ(shift);
}

void SymList::addShift(const Matrix1D<double> &shift)
{
    if (shift.size() != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "SymList::add_shift: Shift vector is not 3x1");
    int i = MAT_YSIZE(__shift);
    __shift.resize(i + 1, 3);
    setShift(i, shift);
}

// Add matrix ==============================================================
void SymList::addMatrices(const Matrix2D<double> &L, const Matrix2D<double> &R,
                          int chain_length)
{
    if (MAT_XSIZE(L) != 4 || MAT_YSIZE(L) != 4 || MAT_XSIZE(R) != 4 || MAT_YSIZE(R) != 4)
        REPORT_ERROR(ERR_MATRIX_SIZE, "SymList::add_matrix: Transformation matrix is not 4x4");
    if (trueSymsNo() == symsNo())
    {
        __L.resize(MAT_YSIZE(__L) + 4, 4);
        __R.resize(MAT_YSIZE(__R) + 4, 4);
        __chain_length.resize(__chain_length.size() + 1);
    }

    setMatrices(true_symNo, L, R);
    __chain_length(__chain_length.size() - 1) = chain_length;
    true_symNo++;
}

// Compute subgroup ========================================================
bool found_not_tried(const Matrix2D<int> &tried, int &i, int &j,
                     int true_symNo)
{
    i = j = 0;
    size_t n = 0;
    while (n != MAT_YSIZE(tried))
    {
        //       if (tried(i, j) == 0 && !(i >= true_symNo && j >= true_symNo))
        if (dMij(tried,i, j) == 0 && !(i >= true_symNo && j >= true_symNo))
            return true;
        if (i != (int)n)
        {
            // Move downwards
            i++;
        }
        else
        {
            // Move leftwards
            j--;
            if (j == -1)
            {
                n++;
                j = n;
                i = 0;
            }
        }
    }
    return false;
}

//#define DEBUG
void SymList::computeSubgroup(double accuracy)
{
    Matrix2D<double> I(4, 4);
    I.initIdentity();
    Matrix2D<double> L1(4, 4), R1(4, 4), L2(4, 4), R2(4, 4), newL(4, 4), newR(4, 4),identity(4,4);
    Matrix2D<int>    tried(true_symNo, true_symNo);
    Matrix1D<double> shift(3);
    shift.initZeros();
    int i, j;
    int new_chain_length;
    identity.initIdentity();
    while (found_not_tried(tried, i, j, true_symNo))
    {
        tried(i, j) = 1;

        // Form new symmetry matrices
        // if (__chain_length(i)+__chain_length(j)>__sym_elements+2) continue;

        getMatrices(i, L1, R1);
        getMatrices(j, L2, R2);
        newL = L1 * L2;
        newR = R1 * R2;
        new_chain_length = __chain_length(i) + __chain_length(j);

        //if (newL.isIdentity() && newR.isIdentity()) continue;
        //rounding error make newR different from identity
        if (newL.equal(identity, accuracy) &&
            newR.equal(identity, accuracy))
            continue;

        // Try to find it in current ones
        bool found;
        found = false;
        for (int l = 0; l < symsNo(); l++)
        {
            getMatrices(l, L1, R1);
            if (newL.equal(L1, accuracy) && newR.equal(R1, accuracy))
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            //#define DEBUG
#ifdef DEBUG
            static int kjhg=0;
            /* std::cout << "Matrix size " << tried.Xdim() << " "
             << "trying " << i << " " << j << " "
             << "chain length=" << new_chain_length << std::endl;
             std::cout << "Result R Sh\n" << newR << shift;
             */
            //std::cerr << "shift" << __shift <<std::endl;
            std::cerr << "newR "  << kjhg++ << "\n" << newR <<std::endl;
#endif

            addMatrices(newL, newR, new_chain_length);
            addShift(shift);
            tried.resize(MAT_YSIZE(tried) + 1, MAT_XSIZE(tried) + 1);
        }
    }
#ifdef DEBUG
    std::cerr << "__R" << __R <<std::endl;
#endif
#undef DEBUG
}
/** Guess Crystallographic space group.
    Return the group number
    http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen. So
    far it has only been implemented for P1 (1), P2_122 (17) ,
    P22_12,
    P4 (75), P4212 (90) and P6 (168) */

int  SymList::crystallographicSpaceGroup(double mag_a, double mag_b,
        double ang_a2b_deg) const
{

    switch (space_group)
    {
    case sym_undefined:
    case sym_P1:
        return(space_group);
    case sym_P4:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P42 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P2_122:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P2_122 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P22_12:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P22_12 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P42_12:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P42_12 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P6:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 120.) > XMIPP_EQUAL_ACCURACY)
        {
            std::cerr << "\nWARNING: marked as P6 but mag_a != mag_b\n"
            << "or ang_a2b !=120" << std::endl;
            std::cerr << "\nWARNING: P1 is assumed\n";
            return(sym_P1);
        }
        else
            return(space_group);
        break;
    default:
        std::cerr << "\n Congratulations: you have found a bug in the\n"
        << "routine crystallographic_space_group or\n"
        << "You have called to this rotuine BEFORE reading\n"
        << "the symmetry info" << std::endl;
        exit(0);
        break;
    }//switch(space_group)  end

}//crystallographicSpaceGroup end
