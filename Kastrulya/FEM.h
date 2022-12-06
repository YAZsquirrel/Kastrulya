#pragma once

//#define DEBUG2
//#define DEBUG

#include "Math_structs.h"
#include <functional>
#include <iomanip>
#include <cmath>
using namespace maths;

class FEM
{
private:
   //real f(knot& knot_) {return 0;}
   //real ug(knot& knot_) {return ;}
   int num_of_knots, num_of_FE, utest;

   real localM[4][4]{};
   real localC[4][4]{};
   real localMb[2][2]{};
   real localG[4][4]{};
   real localGx[4][4]{};
   real localA[4][4]{};

#ifdef DEBUG2
   real f(knot& k, element& e, real time);
   real bound1func(knot &k, real time, int n_mat);
   real bound2func(knot& k, real time, int n_mat);
   real bound3func(knot& k, real time, int n_mat);
   real bound3funcbeta(knot& k, real time, int n_mat);
#endif // DEBUG

   void AddFirstBounds(real time);
   void AddSecondBounds(real time);
   void AddThirdBounds(real time);
   void AddToGlobalMatricies(element& elem);
   void CreateSLAE(bool isTimed);
   void AssembleMatricies(bool isTimed, real time);
   void CreateM(element& elem);
   void CreateG(element& elem);
   void CreateExtraG(element& elem);
   void Createb(element& elem, real time);

   Matrix* A, *M, *G, *Gx;
   std::vector<real> b;
   std::vector<real> d;
   std::vector<real> q1;
   std::vector<real> q2;
   real t_last, dt, tr, u0 = 0;
   Mesh* mesh;


public:
   int GetKnotsNum() { return num_of_knots; }
   int GetHexasNum() { return num_of_FE; }

   std::vector<real>& GetKnots() { return q1; };
   FEM();
   void SolveElliptic();
   void SolveParabolic();
   void Output(std::ofstream& out);

};

