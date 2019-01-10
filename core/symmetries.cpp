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
    int pgGroup = sym_undefined;
    int pgOrder;
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

bool SymList::isSymmetryGroup(FileName fn_sym, int &pgGroup, int &pgOrder)
{
    char G1,G2,G3='\0',G4;
    char auxChar[3];
    //each case check length, check first letter, second, is number
    //Non a point group

    //remove path
    FileName fn_sym_tmp;
    fn_sym_tmp=fn_sym.removeDirectories();
    int mySize=fn_sym_tmp.size();
    bool return_true;
    return_true=false;
    auxChar[2]='\0';
    //size maybe 4 because n maybe a 2 digit number
    if(mySize>4 || mySize<1)
    {
        pgGroup=-1;
        pgOrder=-1;
        return false;
    }
    //get the group character by character
    G1=toupper((fn_sym_tmp.c_str())[0]);
    G2=toupper((fn_sym_tmp.c_str())[1]);
    if (mySize > 2)
    {
        G3=toupper((fn_sym_tmp.c_str())[2]);
        if(mySize > 3)
            G4=toupper((fn_sym.c_str())[3]);
    }
    else
        G4='\0';
    //CN
    if (mySize==2 && G1=='C' && isdigit(G2))
    {
        pgGroup=pg_CN;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    if (mySize==3 && G1=='C' && isdigit(G2) && isdigit(G3))
    {
        pgGroup=pg_CN;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //CI
    else if (mySize==2 && G1=='C' && G2=='I')
    {
        pgGroup=pg_CI;
        pgOrder=-1;
        return_true=true;
    }
    //CS
    else if (mySize==2 && G1=='C' && G2=='S')
    {
        pgGroup=pg_CS;
        pgOrder=-1;
        return_true=true;
    }
    //CNH
    else if (mySize==3 && G1=='C' && isdigit(G2) && G3=='H')
    {
        pgGroup=pg_CNH;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='C' && isdigit(G2) && isdigit(G3) && G4=='H')
    {
        pgGroup=pg_CNH;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //CNV
    else if (mySize==3 && G1=='C' && isdigit(G2) && G3=='V')
    {
        pgGroup=pg_CNV;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='C' && isdigit(G2) && isdigit(G3) && G4=='V')
    {
        pgGroup=pg_CNV;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //SN
    else if (mySize==2 && G1=='S' && isdigit(G2) )
    {
        pgGroup=pg_SN;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==3 && G1=='S' && isdigit(G2) && isdigit(G3) )
    {
        pgGroup=pg_SN;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //DN
    else if (mySize==2 && G1=='D' && isdigit(G2) )
    {
        pgGroup=pg_DN;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    if (mySize==3 && G1=='D' && isdigit(G2) && isdigit(G3))
    {
        pgGroup=pg_DN;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //DNV
    else if (mySize==3 && G1=='D' && isdigit(G2) && G3=='V')
    {
        pgGroup=pg_DNV;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='D' && isdigit(G2) && isdigit(G3) && G4=='V')
    {
        pgGroup=pg_DNV;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //DNH
    else if (mySize==3 && G1=='D' && isdigit(G2) && G3=='H')
    {
        pgGroup=pg_DNH;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='D' && isdigit(G2) && isdigit(G3) && G4=='H')
    {
        pgGroup=pg_DNH;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //T
    else if (mySize==1 && G1=='T')
    {
        pgGroup=pg_T;
        pgOrder=-1;
        return_true=true;
    }
    //TD
    else if (mySize==2 && G1=='T' && G2=='D')
    {
        pgGroup=pg_TD;
        pgOrder=-1;
        return_true=true;
    }
    //TH
    else if (mySize==2 && G1=='T' && G2=='H')
    {
        pgGroup=pg_TH;
        pgOrder=-1;
        return_true=true;
    }
    //O
    else if (mySize==1 && G1=='O')
    {
        pgGroup=pg_O;
        pgOrder=-1;
        return_true=true;
    }
    //OH
    else if (mySize==2 && G1=='O'&& G2=='H')
    {
        pgGroup=pg_OH;
        pgOrder=-1;
        return_true=true;
    }
    //I
    else if (mySize==1 && G1=='I')
    {
        pgGroup=pg_I;
        pgOrder=-1;
        return_true=true;
    }
    //I1
    else if (mySize==2 && G1=='I'&& G2=='1')
    {
        pgGroup=pg_I1;
        pgOrder=-1;
        return_true=true;
    }
    //I2
    else if (mySize==2 && G1=='I'&& G2=='2')
    {
        pgGroup=pg_I2;
        pgOrder=-1;
        return_true=true;
    }
    //I3
    else if (mySize==2 && G1=='I'&& G2=='3')
    {
        pgGroup=pg_I3;
        pgOrder=-1;
        return_true=true;
    }
    //I4
    else if (mySize==2 && G1=='I'&& G2=='4')
    {
        pgGroup=pg_I4;
        pgOrder=-1;
        return_true=true;
    }
    //I5
    else if (mySize==2 && G1=='I'&& G2=='5')
    {
        pgGroup=pg_I5;
        pgOrder=-1;
        return_true=true;
    }
    //IH
    else if (mySize==2 && G1=='I'&& G2=='H')
    {
        pgGroup=pg_IH;
        pgOrder=-1;
        return_true=true;
    }
    //I1H
    else if (mySize==3 && G1=='I'&& G2=='1'&& G3=='H')
    {
        pgGroup=pg_I1H;
        pgOrder=-1;
        return_true=true;
    }
    //I2H
    else if (mySize==3 && G1=='I'&& G2=='2'&& G3=='H')
    {
        pgGroup=pg_I2H;
        pgOrder=-1;
        return_true=true;
    }
    //I3H
    else if (mySize==3 && G1=='I'&& G2=='3'&& G3=='H')
    {
        pgGroup=pg_I3H;
        pgOrder=-1;
        return_true=true;
    }
    //I4H
    else if (mySize==3 && G1=='I'&& G2=='4'&& G3=='H')
    {
        pgGroup=pg_I4H;
        pgOrder=-1;
        return_true=true;
    }
    //I5H
    else if (mySize==3 && G1=='I'&& G2=='5'&& G3=='H')
    {
        pgGroup=pg_I5H;
        pgOrder=-1;
        return_true=true;
    }
    //#define DEBUG7
#ifdef DEBUG7
    std::cerr << "pgGroup" << pgGroup << " pgOrder " << pgOrder << std::endl;
#endif
#undef DEBUG7

    return return_true;
}
void SymList::fillSymmetryClass(const FileName &symmetry, int pgGroup, int pgOrder,
                                std::vector<std::string> &fileContent)
{
    std::ostringstream line1;
    std::ostringstream line2;
    std::ostringstream line3;
    std::ostringstream line4;
    if (pgGroup == pg_CN)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
    }
    else if (pgGroup == pg_CI)
    {
        line1 << "inversion ";
    }
    else if (pgGroup == pg_CS)
    {
        line1 << "mirror_plane 0 0 1";
    }
    else if (pgGroup == pg_CNV)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "mirror_plane 0 1 0";
    }
    else if (pgGroup == pg_CNH)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "mirror_plane 0 0 1";
    }
    else if (pgGroup == pg_SN)
    {
        int order = pgOrder / 2;
        if(2*order != pgOrder)
        {
            std::cerr << "ERROR: order for SN group must be even" << std::endl;
            exit(0);
        }
        line1 << "rot_axis " << order << " 0 0 1";
        line2 << "inversion ";
    }
    else if (pgGroup == pg_DN)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "rot_axis " << "2" << " 1 0 0";
    }
    else if (pgGroup == pg_DNV)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "rot_axis " << "2" << " 1 0 0";
        line3 << "mirror_plane 1 0 0";
    }
    else if (pgGroup == pg_DNH)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "rot_axis " << "2" << " 1 0 0";
        line3 << "mirror_plane 0 0 1";
    }
    else if (pgGroup == pg_T)
    {
        line1 << "rot_axis " << "3" << "  0. 0. 1.";
        line2 << "rot_axis " << "2" << " 0. 0.816496 0.577350";
    }
    else if (pgGroup == pg_TD)
    {
        line1 << "rot_axis " << "3" << "  0. 0. 1.";
        line2 << "rot_axis " << "2" << " 0. 0.816496 0.577350";
        line3 << "mirror_plane 1.4142136 2.4494897 0.0000000";
    }
    else if (pgGroup == pg_TH)
    {
        line1 << "rot_axis " << "3" << "  0. 0. 1.";
        line2 << "rot_axis " << "2" << " 0. -0.816496 -0.577350";
        line3 << "inversion";
    }
    else if (pgGroup == pg_O)
    {
        line1 << "rot_axis " << "3" << "  .5773502  .5773502 .5773502";
        line2 << "rot_axis " << "4" << " 0 0 1";
    }
    else if (pgGroup == pg_OH)
    {
        line1 << "rot_axis " << "3" << "  .5773502  .5773502 .5773502";
        line2 << "rot_axis " << "4" << " 0 0 1";
        line3 << "mirror_plane 0 1 1";
    }
    else if (pgGroup == pg_I || pgGroup == pg_I2)
    {
        line1 << "rot_axis 2  0 0 1";
        line2 << "rot_axis 5  0.525731114  0 0.850650807";
        line3 << "rot_axis 3  0 0.356822076 0.934172364";
    }
    else if (pgGroup == pg_I1)
    {
        line1 << "rot_axis 2  1      0        0";
        line2 << "rot_axis 5 0.85065080702670 0 -0.5257311142635";
        line3 << "rot_axis 3 0.9341723640 0.3568220765 0";
    }
    else if (pgGroup == pg_I3)
    {
        line1 << "rot_axis 2  -0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0. 0. 1.";
        line2 << "rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428";
    }
    else if (pgGroup == pg_I4)
    {
        line1 << "rot_axis 2  0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0.8944271932547096 0 0.4472135909903704";
        line2 << "rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428";
    }
    else if (pgGroup == pg_I5)
    {
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else if (pgGroup == pg_IH || pgGroup == pg_I2H)
    {
        line1 << "rot_axis 2  0 0 1";
        line2 << "rot_axis 5  0.525731114  0 0.850650807";
        line3 << "rot_axis 3  0 0.356822076 0.934172364";
        line4 << "mirror_plane 1 0 0";
    }
    else if (pgGroup == pg_I1H)
    {
        line1 << "rot_axis 2  1      0        0";
        line2 << "rot_axis 5 0.85065080702670 0 -0.5257311142635";
        line3 << "rot_axis 3 0.9341723640 0.3568220765 0";
        line4 << "mirror_plane 0 0 -1";
    }
    else if (pgGroup == pg_I3H)
    {
        line1 << "rot_axis 2  -0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0. 0. 1.";
        line2 << "rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428";
        line4 << "mirror_plane 0.850650807 0  0.525731114";
    }
    else if (pgGroup == pg_I4H)
    {
        line1 << "rot_axis 2  0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0.8944271932547096 0 0.4472135909903704";
        line2 << "rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428";
        line4 << "mirror_plane 0.850650807 0 -0.525731114";
    }
    else if (pgGroup == pg_I5H)
    {
        std::cerr << "ERROR: Symmetry pg_I5H not implemented" << std::endl;
        exit(0);
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }
    if (line1.str().size()>0)
        fileContent.push_back(line1.str());
    if (line2.str().size()>0)
        fileContent.push_back(line2.str());
    if (line3.str().size()>0)
        fileContent.push_back(line3.str());
    if (line4.str().size()>0)
        fileContent.push_back(line4.str());
    //#define DEBUG5
#ifdef DEBUG5

    for (int n=0; n<fileContent.size(); n++)
        std::cerr << fileContent[n] << std::endl;
    std::cerr << "fileContent.size()" << fileContent.size() << std::endl;
#endif
    #undef DEBUG5
}
double SymList::nonRedundantProjectionSphere(int pgGroup, int pgOrder)
{
    if (pgGroup == pg_CN)
    {
        return 4.*PI/pgOrder;
    }
    else if (pgGroup == pg_CI)
    {
        return 4.*PI/2.;
    }
    else if (pgGroup == pg_CS)
    {
        return 4.*PI/2.;
    }
    else if (pgGroup == pg_CNV)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_CNH)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_SN)
    {
        return 4.*PI/pgOrder;
    }
    else if (pgGroup == pg_DN)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_DNV)
    {
        return 4.*PI/pgOrder/4;
    }
    else if (pgGroup == pg_DNH)
    {
        return 4.*PI/pgOrder/4;
    }
    else if (pgGroup == pg_T)
    {
        return 4.*PI/12;
    }
    else if (pgGroup == pg_TD)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_TH)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_O)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_OH)
    {
        return 4.*PI/48;
    }
    else if (pgGroup == pg_I || pgGroup == pg_I2)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I1)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I3)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I4)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I5)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_IH || pgGroup == pg_I2H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I1H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I3H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I4H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I5H)
    {
        return 4.*PI/120;
    }
    else
    {
        std::cerr << "ERROR: Symmetry group, order=" << pgGroup
        << " "
        <<  pgOrder
        << "is not known"
        << std::endl;
        exit(0);
    }
}

void SymList::computeDistance(MetaData &md,
                              bool projdir_mode, bool check_mirrors,
                              bool object_rotation)
{
    MDRow row;
    double rot1, tilt1, psi1;
    double rot2, tilt2, psi2;
    double angDistance;
    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        md.getRow(row,__iter.objId);

        row.getValue(MDL_ANGLE_ROT,rot1);
        row.getValue(MDL_ANGLE_ROT2,rot2);

        row.getValue(MDL_ANGLE_TILT,tilt1);
        row.getValue(MDL_ANGLE_TILT2,tilt2);

        row.getValue(MDL_ANGLE_PSI,psi1);
        row.getValue(MDL_ANGLE_PSI2,psi2);

        angDistance=computeDistance( rot1,  tilt1,  psi1,
                                     rot2,  tilt2,  psi2,
                                     projdir_mode,  check_mirrors,
                                     object_rotation);

        md.setValue(MDL_ANGLE_ROT_DIFF,rot1 - rot2,__iter.objId);
        md.setValue(MDL_ANGLE_TILT_DIFF,tilt1 - tilt2,__iter.objId);
        md.setValue(MDL_ANGLE_PSI_DIFF,psi1 - psi2,__iter.objId);
        md.setValue(MDL_ANGLE_DIFF,angDistance,__iter.objId);
    }

}

double SymList::computeDistance(double rot1, double tilt1, double psi1,
                                double &rot2, double &tilt2, double &psi2,
                                bool projdir_mode, bool check_mirrors,
                                bool object_rotation)
{
    Matrix2D<double> E1, E2;
    Euler_angles2matrix(rot1, tilt1, psi1, E1, false);

    int imax = symsNo() + 1;
    Matrix2D<double>  L(3, 3), R(3, 3);  // A matrix from the list
    double best_ang_dist = 3600;
    double best_rot2=0, best_tilt2=0, best_psi2=0;

    for (int i = 0; i < imax; i++)
    {
        double rot2p, tilt2p, psi2p;
        if (i == 0)
        {
            rot2p = rot2;
            tilt2p = tilt2;
            psi2p = psi2;
        }
        else
        {
            getMatrices(i - 1, L, R, false);
            if (object_rotation)
                Euler_apply_transf(R, L, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
            else
                Euler_apply_transf(L, R, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
        }

        double ang_dist = Euler_distanceBetweenAngleSets_fast(E1,rot2p, tilt2p, psi2p,
                          projdir_mode, E2);

        if (ang_dist < best_ang_dist)
        {
            best_rot2 = rot2p;
            best_tilt2 = tilt2p;
            best_psi2 = psi2p;
            best_ang_dist = ang_dist;
        }

        if (check_mirrors)
        {
        	Euler_mirrorY(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
            double ang_dist_mirror = Euler_distanceBetweenAngleSets_fast(E1,
                                     rot2p, tilt2p, psi2p,projdir_mode, E2);

            if (ang_dist_mirror < best_ang_dist)
            {
                best_rot2 = rot2p;
                best_tilt2 = tilt2p;
                best_psi2 = psi2p;
                best_ang_dist = ang_dist_mirror;
            }

        }
    }
    rot2 = best_rot2;
    tilt2 = best_tilt2;
    psi2 = best_psi2;
    return best_ang_dist;
}

void SymList::breakSymmetry(double rot1, double tilt1, double psi1,
                              double &rot2, double &tilt2, double &psi2
                              )
{
    Matrix2D<double> E1;
    Euler_angles2matrix(rot1, tilt1, psi1, E1, true);
    static bool doRandomize=true;
    Matrix2D<double>  L(3, 3), R(3, 3);  // A matrix from the list

    int i;
    if (doRandomize)
    {
        srand ( time(NULL) );
        doRandomize=false;
    }
    int symOrder = symsNo()+1;
    //std::cerr << "DEBUG_ROB: symOrder: " << symOrder << std::endl;
    i = rand() % symOrder;//59+1
    //std::cerr << "DEBUG_ROB: i: " << i << std::endl;
    if (i < symOrder-1)
    {
        getMatrices(i, L, R);
        //std::cerr  << R << std::endl;
        Euler_matrix2angles(E1 * R, rot2, tilt2, psi2);
    }
    else
    	{
    	//std::cerr << "else" <<std::endl;
    	rot2=rot1; tilt2=tilt1;psi2=psi1;
    	}
//    if (rot2==0)
//:    	std::cerr << "rot2  is zero " << i << R << L << std::endl;
}
