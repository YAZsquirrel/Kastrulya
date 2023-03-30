#pragma once

#include "GridMaker.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <mkl_pardiso.h>

typedef double real;

namespace maths {

   enum MatrixFormat {Dense = 1, SparseRowColumn, SparseProfile, SparseRow};

   struct Matrix {
      std::vector<real> l, u, di, gg;
      std::vector<int> ig, jg;
      std::vector<std::vector<real>> dense;
      size_t dim = 0;
      MatrixFormat format = MatrixFormat::Dense;

      Matrix(MatrixFormat _format = Dense) : format(_format)
      {}
   };

   Matrix* MakeSparseRowColumnFormat(int localsize, int size, Mesh* mesh); // RCF
   Matrix* MakeDenseFormat(int size);
   Matrix* MakeSparseProfileFormat(int localsize, int size, Mesh* mesh);
   void MakeSparseRowFormatFromRCF(Matrix* M, Matrix*& M_srf);
   Matrix* MakeSparseProfileFormatFromRCF(Matrix* M);

   void copy(std::vector<real>& to, std::vector<real>& from);
   void copy(std::vector<int>& to, std::vector<int>& from);
   real scalar(std::vector<real>& v, std::vector<real>& u);
   void AddElement(Matrix* M, int i, int j, real elem);
   void MatxVec(std::vector<real>& v, Matrix* A, std::vector<real>& b);
   void SolveSLAE_LOS(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void SLAEResidualOutput(std::vector<real>& q, maths::Matrix* M, std::vector<real>& b);
   void SolveSLAE_PARDISO(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void SolveSLAE_LOSnKholessky(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void SolveSLAE_Relax(Matrix* M, std::vector<real>& q, std::vector<real>& b, real w);
   void SolveSLAE_LU(Matrix *&LUp, Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void MakeLUFromRCF(Matrix*& LU, maths::Matrix* A);
   void WriteMatrix(Matrix* M);
   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i);
   Matrix* MakeKholessky(Matrix* A);
   void SolveForL(std::vector<real>& q, std::vector<real>& b, Matrix* M);
   void SolveForU(std::vector<real>& q, std::vector<real>& b, Matrix* M);

   void SLAEResidualOutput(std::vector<real>& q, maths::Matrix* M, std::vector<real>& b);

   void ConvertFromSRCFToSRF(Matrix* M_rsf, Matrix* M_csr); 
   void EqualizeRSFToCSR(Matrix* M_srcf, Matrix* M_srf);
}
