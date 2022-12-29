/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include <queue>
#include <fstream>
#include "alglib/linalg.h"
#include "numerical_recipes.h"
#include "matrix2d.h"
#include "bilib/linearalgebra.h"
#include "xmipp_filename.h"
#include "matrix1d.h"
#include "xmipp_funcs.h"
#include <sys/stat.h>

template<typename T>
T Matrix2D<T>::det() const
{
    // (see Numerical Recipes, Chapter 2 Section 5)
    if (mdimx == 0 || mdimy == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "determinant: Matrix is empty");

    if (mdimx != mdimy)
        REPORT_ERROR(ERR_MATRIX_SIZE, "determinant: Matrix is not squared");

    for (size_t i = 0; i < mdimy; i++)
    {
        bool all_zeros = true;
        for (size_t j = 0; j < mdimx; j++)
            if (fabs(MAT_ELEM((*this),i, j)) > XMIPP_EQUAL_ACCURACY)
            {
                all_zeros = false;
                break;
            }

        if (all_zeros)
            return 0;
    }

    // Perform decomposition
    Matrix1D< int > indx;
    T d;
    Matrix2D<T> LU;
    ludcmp(*this, LU, indx, d);

    // Calculate determinant
    for (size_t i = 0; i < mdimx; i++)
        d *= (T) MAT_ELEM(LU,i , i);

    return d;
}

template<typename T>
void Matrix2D<T>::coreInit(const FileName &fn, int Ydim, int Xdim, size_t offset)
{
#ifdef XMIPP_MMAP

    mdimx=Xdim;
    mdimy=Ydim;
    mdim=mdimx*mdimy;
    destroyData=false;
    mappedData=true;
    fdMap = open(fn.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
    if (fdMap == -1)
        REPORT_ERROR(ERR_IO_NOTOPEN,fn);
    const size_t pagesize=sysconf(_SC_PAGESIZE);
    size_t offsetPages=(offset/pagesize)*pagesize;
    size_t offsetDiff=offset-offsetPages;
    if ( (mdataOriginal = (char*) mmap(0,Ydim*Xdim*sizeof(T)+offsetDiff, PROT_READ | PROT_WRITE, MAP_SHARED, fdMap, offsetPages)) == MAP_FAILED )
        REPORT_ERROR(ERR_MMAP_NOTADDR,(String)"mmap failed "+integerToString(errno));
    mdata=(T*)(mdataOriginal+offsetDiff);
#else

    REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

}

template<typename T>
void Matrix2D<T>::read(const FileName &fn)
{
    std::ifstream fhIn;
    fhIn.open(fn.c_str());
    if (!fhIn)
        REPORT_ERROR(ERR_IO_NOTEXIST,fn);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
    fhIn >> MAT_ELEM(*this,i,j);
    fhIn.close();
}

template<typename T>
void Matrix2D<T>::write(const FileName &fn) const
{
    std::ofstream fhOut;
    fhOut.open(fn.c_str());
    if (!fhOut)
        REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"write: Cannot open "+fn+" for output");
    fhOut << *this;
    fhOut.close();
}

#define VIA_BILIB
template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v)
{
    // svdcmp only works with double
    typeCast(a, u);

    // Set size of matrices
    w.initZeros(u.mdimx);
    v.initZeros(u.mdimx, u.mdimx);

    // Call to the numerical recipes routine
#ifdef VIA_NR

    svdcmp(u.mdata,
           u.mdimy, u.mdimx,
           w.vdata,
           v.mdata);
#endif

#ifdef VIA_BILIB

    int status;
    SingularValueDecomposition(u.mdata,
                               u.mdimy, u.mdimx,
                               w.vdata,
                               v.mdata,
                               5000, &status);
#endif
}
#undef VIA_NR
#undef VIA_BILIB

/* Cholesky decomposition -------------------------------------------------- */
void cholesky(const Matrix2D<double> &M, Matrix2D<double> &L)
{
	L=M;
	Matrix1D<double> p;
	p.initZeros(MAT_XSIZE(M));
	choldc(L.adaptForNumericalRecipes2(), MAT_XSIZE(M), p.adaptForNumericalRecipes());
	FOR_ALL_ELEMENTS_IN_MATRIX2D(L)
	if (i==j)
		MAT_ELEM(L,i,j)=VEC_ELEM(p,i);
	else if (i<j)
		MAT_ELEM(L,i,j)=0.0;
}

/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(Matrix2D<double> &u, Matrix1D<double> &w, Matrix2D<double> &v,
            Matrix1D<double> &b, Matrix1D<double> &x)
{
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(u.adaptForNumericalRecipes2(),
           w.adaptForNumericalRecipes(),
           v.adaptForNumericalRecipes2(),
           u.mdimy, u.mdimx,
           b.adaptForNumericalRecipes(),
           x.adaptForNumericalRecipes());
}



void normalizeColumns(Matrix2D<double> &A)
{
	if (MAT_YSIZE(A)<=1)
		return;

	// Compute the mean and standard deviation of each column
	Matrix1D<double> avg, stddev;
	avg.initZeros(MAT_XSIZE(A));
	stddev.initZeros(MAT_XSIZE(A));

	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
	{
		double x=MAT_ELEM(A,i,j);
		VEC_ELEM(avg,j)+=x;
		VEC_ELEM(stddev,j)+=x*x;
	}

	double iN=1.0/MAT_YSIZE(A);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(avg)
	{
        VEC_ELEM(avg,i)*=iN;
        VEC_ELEM(stddev,i)=sqrt(fabs(VEC_ELEM(stddev,i)*iN - VEC_ELEM(avg,i)*VEC_ELEM(avg,i)));
        if (VEC_ELEM(stddev,i)>XMIPP_EQUAL_ACCURACY)
        	VEC_ELEM(stddev,i)=1.0/VEC_ELEM(stddev,i);
        else
        	VEC_ELEM(stddev,i)=0.0;
	}

	// Now normalize
	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
		MAT_ELEM(A,i,j)=(MAT_ELEM(A,i,j)-VEC_ELEM(avg,j))*VEC_ELEM(stddev,j);
}

void normalizeColumnsBetween0and1(Matrix2D<double> &A)
{
	double maxValue,minValue;
	A.computeMaxAndMin(maxValue,minValue);
	double iMaxValue=1.0/(maxValue-minValue);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
		MAT_ELEM(A,i,j)=(MAT_ELEM(A,i,j)-minValue)*iMaxValue;
}

void subtractColumnMeans(Matrix2D<double> &A)
{
	if (MAT_YSIZE(A)<1)
		return;

	// Compute the mean and standard deviation of each column
	Matrix1D<double> avg;
	avg.initZeros(MAT_XSIZE(A));

	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
		VEC_ELEM(avg,j)+=MAT_ELEM(A,i,j);

	double iN=1.0/MAT_YSIZE(A);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(avg)
        VEC_ELEM(avg,i)*=iN;

	// Now normalize
	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
		MAT_ELEM(A,i,j)=MAT_ELEM(A,i,j)-VEC_ELEM(avg,j);
}

void schur(const Matrix2D<double> &M, Matrix2D<double> &O, Matrix2D<double> &T)
{
	alglib::real_2d_array a, s;
	a.setcontent(MAT_YSIZE(M),MAT_XSIZE(M),MATRIX2D_ARRAY(M));
	bool ok=rmatrixschur(a, MAT_YSIZE(M), s);
	if (!ok)
		REPORT_ERROR(ERR_NUMERICAL,"Could not perform Schur decomposition");
	O.resizeNoCopy(M);
	T.resizeNoCopy(M);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
	{
		MAT_ELEM(O,i,j)=s(i,j);
		MAT_ELEM(T,i,j)=a(i,j);
	}
}

void generalizedEigs(const Matrix2D<double> &A, const Matrix2D<double> &B, Matrix1D<double> &D, Matrix2D<double> &P)
{
	int N=(int)MAT_YSIZE(A);
	alglib::real_2d_array a, b, z;
	a.setcontent(N,N,MATRIX2D_ARRAY(A));
	b.setcontent(N,N,MATRIX2D_ARRAY(B));
	alglib::real_1d_array d;
	bool ok=smatrixgevd(a, N, true, b, true, true, 1, d, z);
	if (!ok)
		REPORT_ERROR(ERR_NUMERICAL,"Could not perform eigenvector decomposition");
	D.resizeNoCopy(N);
	memcpy(&VEC_ELEM(D,0),d.getcontent(),N*sizeof(double));
	P.resizeNoCopy(A);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(P)
		MAT_ELEM(P,i,j)=z(i,j);
}

void firstEigs(const Matrix2D<double> &A, size_t M, Matrix1D<double> &D, Matrix2D<double> &P, bool Pneeded)
{
	int N=(int)MAT_YSIZE(A);
	alglib::real_2d_array a, z;
	a.setcontent(N,N,MATRIX2D_ARRAY(A));
	alglib::real_1d_array d;
	bool ok=smatrixevdi(a, N, Pneeded, false, N-M, N-1, d, z);
	if (!ok)
		REPORT_ERROR(ERR_NUMERICAL,"Could not perform eigenvector decomposition");

	D.resizeNoCopy(M);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(D)
		VEC_ELEM(D,i)=d(M-1-i);
	if (Pneeded)
	{
		P.resizeNoCopy(N,M);
		FOR_ALL_ELEMENTS_IN_MATRIX2D(P)
			MAT_ELEM(P,i,j)=z(i,M-1-j);
	}
}

void lastEigs(const Matrix2D<double> &A, size_t M, Matrix1D<double> &D, Matrix2D<double> &P)
{
	int N=(int)MAT_YSIZE(A);
	alglib::real_2d_array a, z;
	a.setcontent(N,N,MATRIX2D_ARRAY(A));
	alglib::real_1d_array d;
	bool ok=smatrixevdi(a, N, true, false, 0, M-1, d, z);
	if (!ok)
		REPORT_ERROR(ERR_NUMERICAL,"Could not perform eigenvector decomposition");

	D.resizeNoCopy(M);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(D)
		VEC_ELEM(D,i)=d(M-1-i);
	memcpy(&VEC_ELEM(D,0),d.getcontent(),M*sizeof(double));
	P.resizeNoCopy(N,M);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(P)
		MAT_ELEM(P,i,j)=z(i,j);
}

void eigsBetween(const Matrix2D<double> &A, size_t I1, size_t I2, Matrix1D<double> &D, Matrix2D<double> &P)
{
	size_t M = I2 - I1 + 1;
	int N=(int)MAT_YSIZE(A);
	alglib::real_2d_array a, z;
	a.setcontent(N,N,MATRIX2D_ARRAY(A));
	alglib::real_1d_array d;

	bool ok=smatrixevdi(a, N, true, false, I1, I2, d, z);
	if (!ok)
		REPORT_ERROR(ERR_NUMERICAL,"Could not perform eigenvector decomposition");

	D.resizeNoCopy(M);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(D)
		VEC_ELEM(D,i)=d(M-1-i);
	memcpy(&VEC_ELEM(D,0),d.getcontent(),M*sizeof(double));
	P.resizeNoCopy(N,M);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(P)
		MAT_ELEM(P,i,j)=z(i,j);
}

void allEigs(const Matrix2D<double> &A, std::vector< std::complex<double> > &eigs)
{
	int N=(int)MAT_YSIZE(A);
	alglib::real_2d_array a, vl, vr;
	a.setcontent(N,N,MATRIX2D_ARRAY(A));
	alglib::real_1d_array wr, wi;

	bool ok=rmatrixevd(a, N, 0, wr, wi, vl, vr);
	if (!ok)
		REPORT_ERROR(ERR_NUMERICAL,"Could not perform eigenvector decomposition");
	eigs.clear();
	for (int n=0; n<N; ++n)
		eigs.push_back(std::complex<double>(wr(n),wi(n)));
}

void connectedComponentsOfUndirectedGraph(const Matrix2D<double> &G, Matrix1D<int> &component)
{
	size_t N=MAT_XSIZE(G);
	component.resizeNoCopy(N);
	component.initConstant(-1);

	int nextComponent=0;
	bool workDone=false;
	std::queue<size_t> toExplore;
	do
	{
		workDone=false;
		// Find next unvisited element
		bool found=false;
		size_t seed=0;
		FOR_ALL_ELEMENTS_IN_MATRIX1D(component)
			if (VEC_ELEM(component,i)<0)
			{
				seed=i;
				found=true;
				break;
			}

		// If found, get its connected component
		if (found)
		{
			int currentComponent=nextComponent;
			nextComponent++;

			VEC_ELEM(component,seed)=currentComponent;
			toExplore.push(seed);
			while (toExplore.size()>0)
			{
				seed=toExplore.front();
				toExplore.pop();
				for (size_t j=seed+1; j<N; ++j)
					if (MAT_ELEM(G,seed,j)>0)
					{
						if (VEC_ELEM(component,j)<0)
						{
							VEC_ELEM(component,j)=currentComponent;
							toExplore.push(j);
						}

					}
			}
			workDone=true;
		}
	} while (workDone);
}

void matrixOperation_AB(const Matrix2D <double> &A, const Matrix2D<double> &B, Matrix2D<double> &C)
{
	C.initZeros(MAT_YSIZE(A), MAT_XSIZE(B));
	for (size_t i = 0; i < MAT_YSIZE(A); ++i)
		for (size_t j = 0; j < MAT_XSIZE(B); ++j)
		{
			double aux=0.;
			for (size_t k = 0; k < MAT_XSIZE(A); ++k)
				aux += MAT_ELEM(A, i, k) * MAT_ELEM(B, k, j);
			MAT_ELEM(C, i, j)=aux;
		}
}

void matrixOperation_Ax(const Matrix2D <double> &A, const Matrix1D<double> &x, Matrix1D<double> &y)
{
    y.initZeros(MAT_YSIZE(A));
	for (size_t i = 0; i < MAT_YSIZE(A); ++i)
	{
		double aux=0.;
		for (size_t k = 0; k < MAT_XSIZE(A); ++k)
			aux += MAT_ELEM(A, i, k) * VEC_ELEM(x, k);
		VEC_ELEM(y, i)=aux;
	}
}

void matrixOperation_AtA(const Matrix2D <double> &A, Matrix2D<double> &B)
{
    B.resizeNoCopy(MAT_XSIZE(A), MAT_XSIZE(A));
    for (size_t i = 0; i < MAT_XSIZE(A); ++i)
        for (size_t j = i; j < MAT_XSIZE(A); ++j)
        {
            double aux=0.;
            for (size_t k = 0; k < MAT_YSIZE(A); ++k)
                aux += MAT_ELEM(A, k, i) * MAT_ELEM(A, k, j);
            MAT_ELEM(B, j, i) = MAT_ELEM(B, i, j) = aux;
        }
}

void matrixOperation_AAt(const Matrix2D <double> &A, Matrix2D<double> &C)
{
	C.initZeros(MAT_YSIZE(A), MAT_YSIZE(A));
	for (size_t i = 0; i < MAT_YSIZE(A); ++i)
		for (size_t j = i; j < MAT_YSIZE(A); ++j)
		{
			double aux=0.;
			for (size_t k = 0; k < MAT_XSIZE(A); ++k)
				aux += MAT_ELEM(A, i, k) * MAT_ELEM(A, j, k);
			MAT_ELEM(C, j, i)=MAT_ELEM(C, i, j)=aux;
		}
}

void matrixOperation_ABt(const Matrix2D <double> &A, const Matrix2D <double> &B, Matrix2D<double> &C)
{
	C.initZeros(MAT_YSIZE(A), MAT_YSIZE(B));
	for (size_t i = 0; i < MAT_YSIZE(A); ++i)
		for (size_t j = 0; j < MAT_YSIZE(B); ++j)
		{
			double aux=0.;
			for (size_t k = 0; k < MAT_XSIZE(A); ++k)
				aux += MAT_ELEM(A, i, k) * MAT_ELEM(B, j, k);
			MAT_ELEM(C, i, j)=aux;
		}
}

void matrixOperation_AtB(const Matrix2D <double> &A, const Matrix2D<double> &B, Matrix2D<double> &C)
{
    C.resizeNoCopy(MAT_XSIZE(A), MAT_XSIZE(B));
    for (size_t i = 0; i < MAT_XSIZE(A); ++i)
        for (size_t j = 0; j < MAT_XSIZE(B); ++j)
        {
            double aux=0.;
            for (size_t k = 0; k < MAT_YSIZE(A); ++k)
                aux += MAT_ELEM(A, k, i) * MAT_ELEM(B, k, j);
			MAT_ELEM(C, i, j)=aux;
        }
}

void matrixOperation_Atx(const Matrix2D <double> &A, const Matrix1D<double> &x, Matrix1D<double> &y)
{
    y.resizeNoCopy(MAT_XSIZE(A));
    for (size_t i = 0; i < MAT_XSIZE(A); ++i)
	{
		double aux=0.;
		for (size_t k = 0; k < MAT_YSIZE(A); ++k)
			aux += MAT_ELEM(A, k, i) * VEC_ELEM(x, k);
		VEC_ELEM(y, i)=aux;
	}
}

void matrixOperation_AtBt(const Matrix2D <double> &A, const Matrix2D<double> &B, Matrix2D<double> &C)
{
	C.initZeros(MAT_XSIZE(A), MAT_YSIZE(B));
	for (size_t i = 0; i < MAT_XSIZE(A); ++i)
		for (size_t j = 0; j < MAT_YSIZE(B); ++j)
		{
			double aux=0.;
			for (size_t k = 0; k < MAT_YSIZE(A); ++k)
				aux += MAT_ELEM(A, k, i) * MAT_ELEM(B, j, k);
			MAT_ELEM(C, i, j)=aux;
		}
}

void matrixOperation_XtAX_symmetric(const Matrix2D<double> &X, const Matrix2D<double> &A, Matrix2D<double> &B)
{
	Matrix2D<double> AX=A*X;
    B.resizeNoCopy(MAT_XSIZE(X), MAT_XSIZE(X));
    for (size_t i = 0; i < MAT_XSIZE(X); ++i)
        for (size_t j = i; j < MAT_XSIZE(X); ++j)
        {
            double aux=0.;
            for (size_t k = 0; k < MAT_YSIZE(X); ++k)
                aux += MAT_ELEM(X, k, i) * MAT_ELEM(AX, k, j);
            MAT_ELEM(B, j, i) = MAT_ELEM(B, i, j) = aux;
        }
}

void matrixOperation_IplusA(Matrix2D<double> &A)
{
	for (size_t i=0; i<MAT_YSIZE(A); ++i)
		MAT_ELEM(A,i,i)+=1;
}

void matrixOperation_IminusA(Matrix2D<double> &A)
{
    FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
        if (i == j)
            MAT_ELEM(A, i, j) = 1 - MAT_ELEM(A, i, j);
        else
            MAT_ELEM(A, i, j) = -MAT_ELEM(A, i, j);
}

void eraseFirstColumn(Matrix2D<double> &A)
{
	Matrix2D<double> Ap;
	Ap.resize(MAT_YSIZE(A),MAT_XSIZE(A)-1);
    for (size_t i = 0; i < MAT_YSIZE(A); ++i)
    	memcpy(&MAT_ELEM(Ap,i,0),&MAT_ELEM(A,i,1),MAT_XSIZE(Ap)*sizeof(double));
    A=Ap;
}

void keepColumns(Matrix2D<double> &A, int j0, int jF)
{
	Matrix2D<double> Ap;
	Ap.resize(MAT_YSIZE(A),jF-j0+1);
    for (size_t i = 0; i < MAT_YSIZE(A); ++i)
    	memcpy(&MAT_ELEM(Ap,i,0),&MAT_ELEM(A,i,j0),MAT_XSIZE(Ap)*sizeof(double));
    A=Ap;
}

void orthogonalizeColumnsGramSchmidt(Matrix2D<double> &M)
{
	for(size_t j1=0; j1<MAT_XSIZE(M); j1++)
	{
		// Normalize column j1
		double norm=0;
		for (size_t i=0; i<MAT_YSIZE(M); i++)
			norm+=MAT_ELEM(M,i,j1)*MAT_ELEM(M,i,j1);
		double K=1.0/sqrt(norm);
		for(size_t i = 0; i<MAT_YSIZE(M); i++)
			MAT_ELEM(M,i,j1) *=K;

		// Update rest of columns
		for(size_t j2=j1+1; j2<MAT_XSIZE(M); j2++)
		{
			// Compute the dot product
			double K = 0;
			for (size_t i=0; i<MAT_YSIZE(M); i++)
				K+=MAT_ELEM(M,i,j1)*MAT_ELEM(M,i,j2);
			for (size_t i=0; i<MAT_YSIZE(M); i++)
				MAT_ELEM(M,i,j2) -= K*MAT_ELEM(M,i,j1);
		}
	}
}

template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    if (VEC_XSIZE(indx)!=A.mdimx)
        indx.resizeNoCopy(A.mdimx);
    ludcmp(LU.adaptForNumericalRecipes2(), A.mdimx,
           indx.adaptForNumericalRecipes(), &d);
}

template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes2(), indx.size(),
           indx.adaptForNumericalRecipes(),
           b.adaptForNumericalRecipes());
}

template<typename T>
void Matrix2D<T>::computeRowMeans(Matrix1D<double> &Xmr) const
{
    Xmr.initZeros(MAT_YSIZE(*this));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        VEC_ELEM(Xmr,i)+=MAT_ELEM(*this,i,j);
    Xmr*=1.0/MAT_XSIZE(*this);
}

template<typename T>
void Matrix2D<T>::computeColMeans(Matrix1D<double> &Xmr) const
{
    Xmr.initZeros(MAT_YSIZE(*this));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        VEC_ELEM(Xmr,j)+=MAT_ELEM(*this,i,j);
    Xmr*=1.0/MAT_YSIZE(*this);
 }

template<>
void Matrix2D<double>::invAlgLib(Matrix2D<double>& result, bool use_lu) const  
{
    if (mdimx < 5) {
        return this->inv(result);
    }

    // copy data to alglib matrix
    const auto N = this->mdimx;
    alglib::real_2d_array a;
    a.setlength(N, N);
    for(size_t i = 0; i < N; ++i) {
        auto dest = a.c_ptr()->ptr.pp_double[i];
        auto src = this->mdata + (i * N);
        memcpy(dest, src, N * sizeof(double));
    }

    // compute inverse
    alglib::ae_int_t info;
    alglib::matinvreport rep;
    if(use_lu) {
        alglib::integer_1d_array pivots;
        alglib::rmatrixlu(a, N, N, pivots);
        alglib::rmatrixluinverse(a, pivots, info, rep);
    } else {
        alglib::rmatrixinverse(a, info, rep); // inplace inverse
    }
    if (1 != (int)info) {
        REPORT_ERROR(ERR_NUMERICAL,"Could not perform matrix inversion using alglib");
    }

    // get result from alglib matrix
    result.resizeNoCopy(*this);
    for(size_t i = 0; i < N; ++i) {
        auto src = a.c_ptr()->ptr.pp_double[i];
        auto dest = result.mdata + (i * N);
        memcpy(dest, src, N * sizeof(double));
    }
}


template<typename T>
void Matrix2D<T>::inv(Matrix2D<T>& result) const
{
    if (mdimx == 0 || mdimy == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "Inverse: Matrix is empty");
    result.initZeros(mdimx, mdimy);
    SPEED_UP_temps0;
    if (mdimx==2)
    {
        M2x2_INV(result,*this);
    }
    else if (mdimx==3)
    {
        M3x3_INV(result,*this);
    }
    else if (mdimx==4)
    {
        M4x4_INV(result,*this);
    }
    else
    {
        // Perform SVD decomposition
        Matrix2D< double > u, v;
        Matrix1D< double > w;
        svdcmp(*this, u, w, v); // *this = U * W * V^t

        double tol = computeMax() * XMIPP_MAX(mdimx, mdimy) * 1e-14;

        // Compute W^-1
        bool invertible = false;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
        {
            if (fabs(VEC_ELEM(w,i)) > tol)
            {
                VEC_ELEM(w,i) = 1.0 / VEC_ELEM(w,i);
                invertible = true;
            }
            else
                VEC_ELEM(w,i) = 0.0;
        }

        if (!invertible)
            return;

        // Compute V*W^-1
        FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
        MAT_ELEM(v,i,j) *= VEC_ELEM(w,j);

        // Compute Inverse
        for (size_t i = 0; i < mdimx; i++)
            for (size_t j = 0; j < mdimy; j++)
                for (size_t k = 0; k < mdimx; k++)
                    MAT_ELEM(result,i,j) += MAT_ELEM(v,i,k) * MAT_ELEM(u,j,k);
    }
}

template<typename T>
void Matrix2D<T>::eigs(Matrix2D<double> &U, Matrix1D<double> &W, Matrix2D<double> &V, Matrix1D<int> &indexes) const
{
  svdcmp(*this, U, W, V);
  indexes.resizeNoCopy(W);
  indexes.enumerate();

  double dAux;
  int iAux;

  FOR_ALL_ELEMENTS_IN_MATRIX1D(W)
  {
      for (int j = i; j > 0 && dMi(W, j) > dMi(W, j-1); --j)
      {
        VEC_SWAP(W, j, j-1, dAux);
        VEC_SWAP(indexes, j, j-1, iAux);
      }
  }
}

template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(const Matrix2D<T>& M)
{
    Matrix1D<T> result;

    if (VEC_XSIZE(*this) != MAT_YSIZE(M))
        REPORT_ERROR(ERR_MATRIX_SIZE, "Not compatible sizes in matrix by vector");

    if (!isRow())
        REPORT_ERROR(ERR_MATRIX_DIM, "Vector is not a row");

    result.initZeros(MAT_XSIZE(M));
    for (size_t j = 0; j < MAT_XSIZE(M); j++)
        for (size_t i = 0; i < MAT_YSIZE(M); i++)
            VEC_ELEM(result,j) += VEC_ELEM(*this,i) * MAT_ELEM(M,i, j);

    result.setRow();
    return result;
}

template<typename T>
Matrix1D<T> Matrix2D<T>::operator*(const Matrix1D<T>& op1) const
{
    Matrix1D<T> result;

    if (mdimx != VEC_XSIZE(op1))
        REPORT_ERROR(ERR_MATRIX_SIZE, "Not compatible sizes in matrix by vector");

    if (!op1.isCol())
        REPORT_ERROR(ERR_MATRIX, "Vector is not a column");

    result.initZeros(mdimy);
    for (size_t i = 0; i < mdimy; i++)
        for (size_t j = 0; j < mdimx; j++)
            VEC_ELEM(result,i) += MAT_ELEM(*this,i, j) * VEC_ELEM(op1,j);

    result.setCol();
    return result;
}

template<typename T>
void Matrix2D<T>::rowSum(Matrix1D<T> &sum) const
{
    sum.initZeros(MAT_YSIZE(*this));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        VEC_ELEM(sum,i)+=MAT_ELEM(*this,i,j);
}

template<typename T>
void Matrix2D<T>::colSum(Matrix1D<T> &sum) const
{
    sum.initZeros(MAT_XSIZE(*this));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        VEC_ELEM(sum,j)+=MAT_ELEM(*this,i,j);
}

template<typename T>
void Matrix2D<T>::rowEnergySum(Matrix1D<T> &sum) const
{
    sum.initZeros(MAT_YSIZE(*this));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        VEC_ELEM(sum,i)+=MAT_ELEM(*this,i,j)*MAT_ELEM(*this,i,j);
}

template<typename T>
void Matrix2D<T>::fromVector(const Matrix1D<T>& op1)
{
    // Null vector => Null matrix
    if (op1.size() == 0)
    {
        clear();
        return;
    }

    // Look at shape and copy values
    if (op1.isRow())
    {
        if (mdimy!=1 || mdimx!=VEC_XSIZE(op1))
            resizeNoCopy(1, VEC_XSIZE(op1));

        for (size_t j = 0; j < VEC_XSIZE(op1); j++)
            MAT_ELEM(*this,0, j) = VEC_ELEM(op1,j);
    }
    else
    {
        if (mdimy!=1 || mdimx!=VEC_XSIZE(op1))
            resizeNoCopy(VEC_XSIZE(op1), 1);

        for (size_t i = 0; i < VEC_XSIZE(op1); i++)
            MAT_ELEM(*this, i, 0) = VEC_ELEM(op1,i);
    }
}

template<typename T>
void Matrix2D<T>::toVector(Matrix1D<T>& op1) const
{
    // Null matrix => Null vector
    if (mdimx == 0 || mdimy == 0)
    {
        op1.clear();
        return;
    }

    // If matrix is not a vector, produce an error
    if (!(mdimx == 1 || mdimy == 1))
        REPORT_ERROR(ERR_MATRIX_DIM,
                     "toVector: Matrix cannot be converted to vector");

    // Look at shape and copy values
    if (mdimy == 1)
    {
        // Row vector
        if (VEC_XSIZE(op1)!=mdimx)
            op1.resizeNoCopy(mdimx);

        memcpy(&VEC_ELEM(op1,0),&MAT_ELEM(*this,0,0),mdimx*sizeof(double));

        op1.setRow();
    }
    else
    {
        // Column vector
        if (VEC_XSIZE(op1)!=mdimy)
            op1.resizeNoCopy(mdimy);

        for (size_t i = 0; i < mdimy; i++)
            VEC_ELEM(op1,i) = MAT_ELEM(*this,i, 0);

        op1.setCol();
    }
}

template<typename T>
void Matrix2D<T>::getRow(size_t i, Matrix1D<T>& v) const
{
    if (mdimx == 0 || mdimy == 0)
    {
        v.clear();
        return;
    }

    if (i >= mdimy)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "getRow: Matrix subscript (i) greater than matrix dimension");

    if (VEC_XSIZE(v)!=mdimx)
        v.resizeNoCopy(mdimx);
    memcpy(&VEC_ELEM(v,0),&MAT_ELEM(*this,i,0),mdimx*sizeof(T));

    v.setRow();
}

template<typename T>
void Matrix2D<T>::getCol(size_t j, Matrix1D<T>& v) const
{
    if (mdimx == 0 || mdimy == 0)
    {
        v.clear();
        return;
    }

    if (j >= mdimx)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"getCol: Matrix subscript (j) greater than matrix dimension");

    if (VEC_XSIZE(v)!=mdimy)
        v.resizeNoCopy(mdimy);
    for (size_t i = 0; i < mdimy; i++)
        VEC_ELEM(v,i) = MAT_ELEM(*this,i, j);

    v.setCol();
}

template<typename T>
void Matrix2D<T>::setRow(size_t i, const Matrix1D<T>& v)
{
    if (mdimx == 0 || mdimy == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "setRow: Target matrix is empty");

    if (i >= mdimy)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setRow: Matrix subscript (i) out of range");

    if (VEC_XSIZE(v) != mdimx)
        REPORT_ERROR(ERR_MATRIX_SIZE,
                     "setRow: Vector dimension different from matrix one");

    if (!v.isRow())
        REPORT_ERROR(ERR_MATRIX_DIM, "setRow: Not a row vector in assignment");

    memcpy(&MAT_ELEM(*this,i,0),&VEC_ELEM(v,0),mdimx*sizeof(double));
}

template<typename T>
void Matrix2D<T>::setCol(size_t j, const Matrix1D<T>& v)
{
    if (mdimx == 0 || mdimy == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "setCol: Target matrix is empty");

    if (j>= mdimx)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setCol: Matrix subscript (j) out of range");

    if (VEC_XSIZE(v) != mdimy)
        REPORT_ERROR(ERR_MATRIX_SIZE,
                     "setCol: Vector dimension different from matrix one");

    if (!v.isCol())
        REPORT_ERROR(ERR_MATRIX_DIM, "setCol: Not a column vector in assignment");

    for (size_t i = 0; i < mdimy; i++)
        MAT_ELEM(*this,i, j) = VEC_ELEM(v,i);
}

template<typename T>
void Matrix2D<T>::getDiagonal(Matrix1D<T> &d) const
{
    d.resizeNoCopy(MAT_XSIZE(*this));
    for (size_t i=0; i<MAT_XSIZE(*this); ++i)
        VEC_ELEM(d,i)=MAT_ELEM(*this,i,i);
}

template<typename T>
Matrix2D<T>& Matrix2D<T>::operator=(const Matrix2D<T>& op1)
{
   if (&op1 != this)
   {
       if (MAT_XSIZE(*this)!=MAT_XSIZE(op1) ||
           MAT_YSIZE(*this)!=MAT_YSIZE(op1))
           resizeNoCopy(op1);
       memcpy(mdata,op1.mdata,op1.mdim*sizeof(T));
   }

   return *this;
}

template<typename T>
void Matrix2D<T>::coreInit()
{
   mdimx=mdimy=mdim=0;
   mdata=NULL;
   mdataOriginal=NULL;
   destroyData=true;
   mappedData=false;
   fdMap=-1;
}

template<typename T>
void Matrix2D<T>::coreAllocate( int _mdimy, int _mdimx)
{
   if (_mdimy <= 0 ||_mdimx<=0)
   {
       clear();
       return;
   }

   mdimx=_mdimx;
   mdimy=_mdimy;
   mdim=_mdimx*_mdimy;
   mdata = new T [mdim];
   mdataOriginal = NULL;
   mappedData=false;
   fdMap=-1;
   if (mdata == NULL)
       REPORT_ERROR(ERR_MEM_NOTENOUGH, "coreAllocate: No space left");
}

template<typename T>
void Matrix2D<T>::coreDeallocate()
{
   if (mdata != NULL && destroyData)
       delete[] mdata;
   if (mappedData)
   {
#ifdef XMIPP_MMAP
       munmap(mdataOriginal,mdimx*mdimy*sizeof(T));
       close(fdMap);
#else

       REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

   }
   mdata=NULL;
   mdataOriginal=NULL;
}

template<typename T>
void Matrix2D<T>::resize(size_t Ydim, size_t Xdim, bool noCopy)
{

   if (Xdim == mdimx && Ydim == mdimy)
       return;

   if (Xdim <= 0 || Ydim <= 0)
   {
       clear();
       return;
   }

   T * new_mdata;
   size_t YXdim=Ydim*Xdim;

   try
   {
       new_mdata = new T [YXdim];
   }
   catch (std::bad_alloc &)
   {
       REPORT_ERROR(ERR_MEM_NOTENOUGH, "Allocate: No space left");
   }

   // Copy needed elements, fill with 0 if necessary
   if (!noCopy)
   {
       T zero=0; // Useful for complexes
       for (size_t i = 0; i < Ydim; i++)
           for (size_t j = 0; j < Xdim; j++)
           {
               T *val=NULL;
               if (i >= mdimy)
                   val = &zero;
               else if (j >= mdimx)
                   val = &zero;
               else
                   val = &mdata[i*mdimx + j];
               new_mdata[i*Xdim+j] = *val;
           }
   }
   else
       memset(new_mdata,0,YXdim*sizeof(T));

   // deallocate old vector
   coreDeallocate();

   // assign *this vector to the newly created
   mdata = new_mdata;
   mdimx = Xdim;
   mdimy = Ydim;
   mdim = Xdim * Ydim;
   mappedData = false;
}


template<typename T>
void Matrix2D<T>::mapToFile(const FileName &fn, int Ydim, int Xdim, size_t offset)
{
   if (mdata!=NULL)
       clear();

#ifdef XMIPP_MMAP

   coreInit(fn,Ydim,Xdim,offset);
#else

   resizeNoCopy(Ydim, Xdim);
#endif

}

template<typename T>
void Matrix2D<T>::submatrix(int i0, int j0, int iF, int jF)
{
   if (i0 < 0 || j0 < 0 || iF >= MAT_YSIZE(*this) || jF >= MAT_XSIZE(*this))
       REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"Submatrix indexes out of bounds");
   Matrix2D<T> result(iF - i0 + 1, jF - j0 + 1);

   FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
   MAT_ELEM(result, i, j) = MAT_ELEM(*this, i+i0, j+j0);

   *this = result;
}

template<typename T>
void Matrix2D<T>::initRandom(size_t Ydim, size_t Xdim, double op1, double op2, RandomMode mode)
{
   if (mdimx!=Xdim || mdimy!=Ydim)
       resizeNoCopy(Ydim, Xdim);
   for (size_t j = 0; j < mdim; j++)
      mdata[j] = static_cast< T > (mode == RND_UNIFORM ? rnd_unif(op1, op2) : rnd_gaus(op1, op2));
}

/** Initialize to gaussian numbers */
template<typename T>
void Matrix2D<T>::initGaussian(int Ydim, int Xdim, double op1, double op2)
{
 initRandom(Ydim, Xdim, op1, op2, RND_GAUSSIAN);
}


template<typename T>
void Matrix2D<T>::initGaussian(int dim, double var)
{
   double center = ((double)dim)/2;
   initZeros(dim, dim);
   for (int i = 0; i < dim; i++)
       for (int j = 0; j < dim; j++)
           MAT_ELEM(*this,i,j) = std::exp(-( (i-center)*(i-center)+(j-center)*(j-center) )/(2*var*var));
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator*(const Matrix2D<T>& op1) const
{
   Matrix2D<T> result;
   if (mdimx != op1.mdimy)
       REPORT_ERROR(ERR_MATRIX_SIZE, "Not compatible sizes in matrix multiplication");

   result.initZeros(mdimy, op1.mdimx);
   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < op1.mdimx; j++)
           for (size_t k = 0; k < mdimx; k++)
               MAT_ELEM(result,i, j) += MAT_ELEM(*this,i, k) * MAT_ELEM(op1, k, j);
   return result;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator+(const Matrix2D<T>& op1) const
{
   Matrix2D<T> result;
   if (mdimx != op1.mdimx || mdimy != op1.mdimy)
       REPORT_ERROR(ERR_MATRIX_SIZE, "operator+: Not same sizes in matrix summation");

   result.initZeros(mdimy, mdimx);
   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           result(i, j) = (*this)(i, j) + op1(i, j);

   return result;
}

template<typename T>
void Matrix2D<T>::operator+=(const Matrix2D<T>& op1) const
{
   if (mdimx != op1.mdimx || mdimy != op1.mdimy)
       REPORT_ERROR(ERR_MATRIX_SIZE, "operator+=: Not same sizes in matrix summation");

   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           MAT_ELEM(*this,i, j) += MAT_ELEM(op1, i, j);
}

template<typename T>
Matrix2D<T> Matrix2D<T>::operator-(const Matrix2D<T>& op1) const
{
   Matrix2D<T> result;
   if (mdimx != op1.mdimx || mdimy != op1.mdimy)
       REPORT_ERROR(ERR_MATRIX_SIZE, "operator-: Not same sizes in matrix summation");

   result.initZeros(mdimy, mdimx);
   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           result(i, j) = (*this)(i, j) - op1(i, j);

   return result;
}

template<typename T>
void Matrix2D<T>::operator-=(const Matrix2D<T>& op1) const
{
   if (mdimx != op1.mdimx || mdimy != op1.mdimy)
       REPORT_ERROR(ERR_MATRIX_SIZE, "operator-=: Not same sizes in matrix summation");

   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           MAT_ELEM(*this,i, j) -= MAT_ELEM(op1, i, j);
}

template<typename T>
bool Matrix2D<T>::equal(const Matrix2D<T>& op,
          double accuracy) const
{
   if (!sameShape(op))
       return false;
   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           if (fabs( MAT_ELEM(*this,i,j) - MAT_ELEM(op,i,j) ) > accuracy)
           {
               //std::cerr << "DEBUG_ROB: MAT_ELEM(*this,i,j): " << MAT_ELEM(*this,i,j) << std::endl;
               //std::cerr << "DEBUG_ROB: MAT_ELEM(op,i,j): " << MAT_ELEM(op,i,j) << std::endl;
               return false;
           }
   return true;
}

template<typename T>
bool Matrix2D<T>::equalAbs(const Matrix2D<T>& op,
          double accuracy) const
{
   if (!sameShape(op))
       return false;
   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           if ( (fabs( MAT_ELEM(*this,i,j)) - fabs(MAT_ELEM(op,i,j)) )> accuracy)
           {
               //std::cerr << "DEBUG_ROB: MAT_ELEM(*this,i,j): " << MAT_ELEM(*this,i,j) << std::endl;
               //std::cerr << "DEBUG_ROB: MAT_ELEM(op,i,j): " << MAT_ELEM(op,i,j) << std::endl;
               return false;
           }
   return true;
}

template<typename T>
T Matrix2D<T>::computeMax() const
{
   if (mdim <= 0)
       return static_cast< T >(0);

   T maxval = mdata[0];
   for (size_t n = 0; n < mdim; n++)
       if (mdata[n] > maxval)
           maxval = mdata[n];
   return maxval;
}

template<typename T>
T Matrix2D<T>::computeMin() const
{
   if (mdim <= 0)
       return static_cast< T >(0);

   T minval = mdata[0];
   for (size_t n = 0; n < mdim; n++)
       if (mdata[n] < minval)
           minval = mdata[n];
   return minval;
}

template<typename T>
void Matrix2D<T>::computeMaxAndMin(T &maxValue, T &minValue) const
{
   maxValue=minValue=0;
   if (mdim <= 0)
       return;

   maxValue = minValue = mdata[0];
   for (size_t n = 0; n < mdim; n++)
   {
       T val=mdata[n];
       if (val < minValue)
           minValue = val;
       else if (val > maxValue)
           maxValue = val;
   }
}

template<typename T>
T** Matrix2D<T>::adaptForNumericalRecipes() const
{
   T** m = NULL;
   ask_Tmatrix(m, 1, mdimy, 1, mdimx);

   for (int i = 0; i < mdimy; i++)
       for (int j = 0; j < mdimx; j++)
           m[i+1][j+1] = mdata[i*mdimx + j];

   return m;
}

template<typename T>
void Matrix2D<T>::loadFromNumericalRecipes(T** m, int Ydim, int Xdim)
{
   if (mdimx!=Xdim || mdimy!=Ydim)
       resizeNoCopy(Ydim, Xdim);

   for (int i = 1; i <= Ydim; i++)
       for (int j = 1; j <= Xdim; j++)
           MAT_ELEM(*this,i - 1, j - 1) = m[i][j];
}

template<typename T>
T Matrix2D<T>::trace() const
{
   size_t d=std::min(MAT_XSIZE(*this),MAT_YSIZE(*this));
   T retval=0;
   for (size_t i=0; i<d; ++i)
       retval+=MAT_ELEM(*this,i,i);
   return retval;
}

template<typename T>
Matrix2D<T> Matrix2D<T>::transpose() const
{
   Matrix2D<T> result(mdimx, mdimy);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
   MAT_ELEM(result,i,j) = MAT_ELEM((*this),j,i);
   return result;
}

template<typename T>
bool Matrix2D<T>::isIdentity() const
{
   for (size_t i = 0; i < mdimy; i++)
       for (size_t j = 0; j < mdimx; j++)
           if (i != j)
           {
               if (MAT_ELEM(*this,i,j)!=0)
                   return false;
           }
           else
           {
               if (MAT_ELEM(*this,i,j)!=1)
                   return false;
           }
   return true;
}

template void ludcmp<double>(Matrix2D<double> const&, Matrix2D<double>&, Matrix1D<int>&, double&);
template void lubksb<double>(Matrix2D<double> const&, Matrix1D<int>&, Matrix1D<double>&);
template void svdcmp<float>(Matrix2D<float> const&, Matrix2D<double>&, Matrix1D<double>&, Matrix2D<double>&);
template void svdcmp<double>(Matrix2D<double> const&, Matrix2D<double>&, Matrix1D<double>&, Matrix2D<double>&);
template Matrix1D<double> Matrix1D<double>::operator*(Matrix2D<double> const&);
template class Matrix2D<double>;
template class Matrix2D<int>;
template class Matrix2D<float>;
template class Matrix2D<unsigned char>;
