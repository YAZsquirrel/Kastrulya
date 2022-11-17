#pragma once

#include "GridMaker.h"
#include <iostream>
#include <fstream>
#include <vector>
typedef double real;

namespace maths {

   struct Matrix {
      std::vector<real> l, u, di;
      std::vector<int> ig, jg;
      size_t dim = 0;
   };



   Matrix* MakeSparseFormat(int localsize, int elemsize, Mesh* mesh);
   void copy(std::vector<real>& to, std::vector<real>& from);
   void copy(std::vector<int>& to, std::vector<int>& from);
   real scalar(std::vector<real>& v, std::vector<real>& u);
   void AddElement(Matrix* M, int i, int j, real elem);
   void MatxVec(std::vector<real>& v, Matrix* A, std::vector<real>& b);
   void SolveSLAE(Matrix* M, std::vector<real>& q, std::vector<real>& b);
   void WriteMatrix(Matrix* M);
   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i);
}