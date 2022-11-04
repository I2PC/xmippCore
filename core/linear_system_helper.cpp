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

#include "linear_system_helper.h"
#include <random>
#include <algorithm>
#include <cassert>

// Solve linear systems ---------------------------------------------------
 // FIXME deprecated (use solveLinearSystem(WeightedLeastSquaresHelperMany &h, std::vector<Matrix1D<double>> &results))
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

void solveLinearSystem(WeightedLeastSquaresHelperMany &h, std::vector<Matrix1D<double>> &results)
{
    const size_t sizeI = h.A.mdimy;
    const size_t sizeJ = h.A.mdimx;
    h.AtA.resizeNoCopy(sizeJ, sizeJ);
    h.Atb.resizeNoCopy(sizeJ);
    // compute AtA
    h.At = h.A.transpose();
    // for each row of the output
    for (size_t i = 0; i < sizeJ; ++i)
    {
        // for each column of the output (notice that we read row-wise from At)
        for (size_t j = i; j < sizeJ; ++j)
        {
            double AtA_ij = 0;
            // multiply elementwise row by row from At
            for (size_t k = 0; k < sizeI; ++k) {
                AtA_ij += h.At.mdata[i * sizeI + k] * h.At.mdata[j * sizeI + k];
            }
            // AtA is symmetric, so we don't need to iterate over everything
            h.AtA.mdata[i * sizeJ + j] = h.AtA.mdata[j * sizeJ + i] = AtA_ij;
        }
    }
    // Compute the inverse of AtA
    h.AtA.inv(h.AtAinv); // AtAinv is also square matrix, same size as AtA

    assert(results.size() == h.bs.size());
    auto res_it = results.begin();
    auto b_it = h.bs.begin(); 
    for (; res_it != results.end(); ++res_it, ++b_it) {
        for (size_t i = 0; i < sizeJ; ++i)
        {
            double Atb_i = 0;
            for (size_t k = 0; k < sizeI; ++k) {
                Atb_i += h.A.mdata[k * sizeJ + i] * b_it->vdata[k];
            }
            h.Atb.vdata[i] = Atb_i;
        }
        // Now multiply by Atb
        res_it->initZeros(sizeJ);
        for (size_t i = 0; i < sizeJ; i++) { 
            for (size_t j = 0; j < sizeJ; j++) {
                res_it->vdata[i] += h.AtAinv.mdata[i * sizeJ +j] * h.Atb.vdata[j];
            }
        }
    }
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

void weightedLeastSquares(WeightedLeastSquaresHelperMany &h, std::vector<Matrix1D<double>> &results)
{
    const size_t sizeW = h.w.vdim;
    h.w_sqrt.resizeNoCopy(h.w);
    // compute weights
    for(size_t i = 0; i < sizeW; ++i) {
        h.w_sqrt.vdata[i] = sqrt(h.w.vdata[i]);
    }
    // update the matrix
    const size_t sizeX = h.A.mdimx;
    for(size_t i = 0; i < sizeW; ++i) {
        const size_t offset = i * sizeX; 
        const auto v = h.w_sqrt[i];
        for (size_t j = 0; j < sizeX; ++j)
            h.A.mdata[offset + j] *= v;
    }
    // update values
    for (auto &b : h.bs) {
        for(size_t i = 0; i < sizeW; ++i) {
            b.vdata[i] *= h.w_sqrt[i];
        }
    }
    solveLinearSystem(h, results);
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
