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

#include <algorithm>
#include <queue>
#include <random>
#include <fstream>
#include "alglib/ap.h"
#include "alglib/linalg.h"
#include "numerical_recipes.h"
#include "matrix2d.h"
#include "bilib/linearalgebra.h"
#include "xmipp_filename.h"

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

// Solve linear systems ---------------------------------------------------
void solveLinearSystem(PseudoInverseHelper &h, Matrix1D<double> &result)
{
	Matrix2D<double> &A=h.A;
	Matrix1D<double> &b=h.b;
	Matrix2D<double> &AtA=h.AtA;
	Matrix2D<double> &AtAinv=h.AtAinv;
	Matrix1D<double> &Atb=h.Atb;

	// Compute AtA and Atb
	int I=MAT_YSIZE(A);
	int J=MAT_XSIZE(A);
	AtA.initZeros(J,J);
	Atb.initZeros(J);
	for (int i=0; i<J; ++i)
	{
		for (int j=0; j<J; ++j)
		{
			double AtA_ij=0;
			for (int k=0; k<I; ++k)
				AtA_ij+=MAT_ELEM(A,k,i)*MAT_ELEM(A,k,j);
			MAT_ELEM(AtA,i,j)=AtA_ij;
		}
		double Atb_i=0;
		for (int k=0; k<I; ++k)
			Atb_i+=MAT_ELEM(A,k,i)*VEC_ELEM(b,k);
		VEC_ELEM(Atb,i)=Atb_i;
	}

	// Compute the inverse of AtA
	AtA.inv(AtAinv);

	// Now multiply by Atb
	result.initZeros(J);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(AtAinv)
		VEC_ELEM(result,i)+=MAT_ELEM(AtAinv,i,j)*VEC_ELEM(Atb,j);
}

// Solve linear systems ---------------------------------------------------
void weightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result)
{
	Matrix2D<double> &A=h.A;
	Matrix1D<double> &b=h.b;
	Matrix1D<double> &w=h.w;

	// See http://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares
	FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
	{
		double wii=sqrt(VEC_ELEM(w,i));
		VEC_ELEM(b,i)*=wii;
		for (size_t j=0; j<MAT_XSIZE(A); ++j)
			MAT_ELEM(A,i,j)*=wii;
	}
	solveLinearSystem(h,result);
}

// Solve linear system with RANSAC ----------------------------------------
//#define DEBUG
//#define DEBUG_MORE
double ransacWeightedLeastSquaresBasic(WeightedLeastSquaresHelper &h, Matrix1D<double> &result,
		double tol, int Niter, double outlierFraction)
{
	int N=MAT_YSIZE(h.A); // Number of equations
	int M=MAT_XSIZE(h.A); // Number of unknowns

	// Initialize a vector with all equation indexes
	std::vector<int> eqIdx;
	eqIdx.reserve(N);
	for (int n=0; n<N; ++n)
		eqIdx.push_back(n);
	int *eqIdxPtr=&eqIdx[0];

#ifdef DEBUG_MORE
	// Show all equations
	for (int n=0; n<N; n++)
	{
		std::cout << "Eq. " << n << " w=" << VEC_ELEM(h.w,n) << " b=" << VEC_ELEM(h.b,n) << " a=";
		for (int j=0; j<M; j++)
			std::cout << MAT_ELEM(h.A,n,j) << " ";
		std::cout << std::endl;
	}
#endif

	// Resize a WLS helper for solving the MxM equation systems
	PseudoInverseHelper haux;
	haux.A.resizeNoCopy(M,M);
	haux.b.resizeNoCopy(M);
	Matrix2D<double> &A=haux.A;
	Matrix1D<double> &b=haux.b;

	// Solve Niter randomly chosen equation systems
	double bestError=1e38;
	const int Mdouble=M*sizeof(double);
	int minNewM=(int)((1.0-outlierFraction)*N-M);
	if (minNewM<0)
		minNewM=0;

	Matrix1D<double> resultAux;
	Matrix1D<int> idxIn(N);
	WeightedLeastSquaresHelper haux2;
	for (int it=0; it<Niter; ++it)
	{
#ifdef DEBUG_MORE
		std::cout << "Iter: " << it << std::endl;
#endif
        idxIn.initZeros();

        // Randomly select M equations
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(eqIdx.begin(), eqIdx.end(), g);

        // Select the equation system
        for (int i=0; i<M; ++i)
        {
        	int idx=eqIdxPtr[i];
        	memcpy(&MAT_ELEM(A,i,0),&MAT_ELEM(h.A,idx,0),Mdouble);
        	VEC_ELEM(b,i)=VEC_ELEM(h.b,idx);
        	VEC_ELEM(idxIn,idx)=1;
#ifdef DEBUG_MORE
		std::cout << "    Using Eq.: " << idx << " for first solution" << std::endl;
#endif
        }

        // Solve the equation system
        // We use LS because the weight of some of the equations might be too low
        // and then the system is ill conditioned
        solveLinearSystem(haux, resultAux);

        // Study the residuals of the rest
        int newM=0;
        for (int i=M+1; i<N; ++i)
        {
        	int idx=eqIdxPtr[i];
        	double bp=0;
        	for (int j=0; j<M; ++j)
        		bp+=MAT_ELEM(h.A,idx,j)*VEC_ELEM(resultAux,j);
        	if (fabs(bp-VEC_ELEM(h.b,idx))<tol)
        	{
        		VEC_ELEM(idxIn,idx)=1;
        		++newM;
        	}
#ifdef DEBUG_MORE
		std::cout << "    Checking Eq.: " << idx << " err=" << bp-VEC_ELEM(h.b,idx) << std::endl;
#endif
        }

        // If the model represent more points
        if (newM>minNewM)
        {
        	Matrix2D<double> &A2=haux2.A;
        	Matrix1D<double> &b2=haux2.b;
        	Matrix1D<double> &w2=haux2.w;
        	A2.resizeNoCopy(M+newM,M);
        	b2.resizeNoCopy(M+newM);
        	w2.resizeNoCopy(M+newM);

            // Select the equation system
        	int targeti=0;
            for (int i=0; i<N; ++i)
            	if (VEC_ELEM(idxIn,i))
            	{
					memcpy(&MAT_ELEM(A2,targeti,0),&MAT_ELEM(h.A,i,0),Mdouble);
					VEC_ELEM(b2,targeti)=VEC_ELEM(h.b,i);
					VEC_ELEM(w2,targeti)=VEC_ELEM(h.w,i);
					++targeti;
            	}

            // Solve it with WLS
            weightedLeastSquares(haux2, resultAux);

            // Compute the mean error
            double err=0;
            for (int i=0; i<M+newM; ++i)
			{
				double bp=0;
				for (int j=0; j<M; ++j)
					bp+=MAT_ELEM(A2,i,j)*VEC_ELEM(resultAux,j);
				err+=fabs(VEC_ELEM(b2,i)-bp)*VEC_ELEM(w2,i);
			}
            err/=(M+newM);
            if (err<bestError)
            {
            	bestError=err;
            	result=resultAux;
#ifdef DEBUG
            	std::cout << "Best solution iter: " << it << " Error=" << err << " frac=" << (float)(M+newM)/VEC_XSIZE(h.b) << std::endl;
#ifdef DEBUG_MORE
            	std::cout << "Result:" << result << std::endl;
                for (int i=0; i<M+newM; ++i)
    			{
    				double bp=0;
    				for (int j=0; j<M; ++j)
    					bp+=MAT_ELEM(A2,i,j)*VEC_ELEM(resultAux,j);
    				std::cout << "Eq. " << i << " w=" << VEC_ELEM(w2,i) << " b2=" << VEC_ELEM(b2,i) << " bp=" << bp << std::endl;
    				err+=fabs(VEC_ELEM(b2,i)-bp)*VEC_ELEM(w2,i);
    			}
#endif
#endif
            }
        }
	}
	return bestError;
}
#undef DEBUG

struct ThreadRansacArgs {
	// Input
	int myThreadID;
	WeightedLeastSquaresHelper * h;
	double tol;
	int Niter;
	double outlierFraction;

	// Output
	Matrix1D<double> result;
	double error;
};

void * threadRansacWeightedLeastSquares(void * args)
{
	ThreadRansacArgs * master = (ThreadRansacArgs *) args;
	master->error=ransacWeightedLeastSquaresBasic(*(master->h), master->result,
			master->tol, master->Niter, master->outlierFraction);
	return NULL;
}

void ransacWeightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result,
		double tol, int Niter, double outlierFraction, int Nthreads)
{
	// Read and preprocess the images
	pthread_t * th_ids = new pthread_t[Nthreads];
	ThreadRansacArgs * th_args = new ThreadRansacArgs[Nthreads];
	for (int nt = 0; nt < Nthreads; nt++) {
		// Passing parameters to each thread
		th_args[nt].myThreadID = nt;
		th_args[nt].h = &h;
		th_args[nt].tol = tol;
		th_args[nt].Niter = Niter/Nthreads;
		th_args[nt].outlierFraction = outlierFraction;
		pthread_create((th_ids + nt), NULL, threadRansacWeightedLeastSquares,
				(void *) (th_args + nt));
	}

	// Waiting for threads to finish
	double err=1e38;
	for (int nt = 0; nt < Nthreads; nt++)
	{
		pthread_join(*(th_ids + nt), NULL);
		if (th_args[nt].error<err)
		{
			err=th_args[nt].error;
			result=th_args[nt].result;
		}
	}

    // Threads structures are not needed any more
    delete []th_ids;
    delete []th_args;
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

template void ludcmp<double>(Matrix2D<double> const&, Matrix2D<double>&, Matrix1D<int>&, double&);
template void lubksb<double>(Matrix2D<double> const&, Matrix1D<int>&, Matrix1D<double>&);
template void svdcmp<float>(Matrix2D<float> const&, Matrix2D<double>&, Matrix1D<double>&, Matrix2D<double>&);
template void Matrix2D<int>::write(FileName const&) const;
template class Matrix2D<double>;
