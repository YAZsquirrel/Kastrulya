#pragma once

#include "GridMaker.h"
#include <iostream>
#include <fstream>
#include <vector>

typedef double real;

namespace maths {

   enum MatrixFormat {Dense = 1, SparseRowColumn, SparseProfile};
   struct Matrix {
      std::vector<real> l, u, di;
      std::vector<int> ig, jg;
      std::vector<std::vector<real>> dense;
      size_t dim = 0;
      MatrixFormat format = MatrixFormat::Dense;
   };

   

   Matrix* MakeSparseRowColumnFormat(int localsize, int size, Mesh* mesh); // RCF
   Matrix* MakeDenseFormat(int size);
   Matrix* MakeSparseProfileFormat(int localsize, int size, Mesh* mesh);
   Matrix* MakeSparseProfileFormatFromRCF(Matrix* M);

   void copy(std::vector<real>& to, std::vector<real>& from);
   void copy(std::vector<int>& to, std::vector<int>& from);
   real scalar(std::vector<real>& v, std::vector<real>& u);
   void AddElement(Matrix* M, int i, int j, real elem);
   void MatxVec(std::vector<real>& v, Matrix* A, std::vector<real>& b);
   void SolveSLAE_LOS(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void SolveSLAE_LOSnKholessky(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void SolveSLAE_Relax(Matrix* M, std::vector<real>& q, std::vector<real>& b, real w);
   void SolveSLAE_LU(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void WriteMatrix(Matrix* M);
   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i);
   Matrix* MakeKholessky(Matrix* A);
   void SolveForL(std::vector<real>& q, std::vector<real>& b, Matrix* M);
   void SolveForU(std::vector<real>& q, std::vector<real>& b, Matrix* M);
}