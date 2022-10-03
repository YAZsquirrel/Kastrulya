#pragma once
#include "Math_structs.h"
//#include <functional>
#include <iomanip>
#include <cmath>
using namespace maths;

class FEM
{
private:
   real f(knot& knot_) {return 0;}
   real ug(knot& knot_) {return u0;}
   int num_of_knots, num_of_FE, un;

   real localM[4][4]; 
   real localMb[2][2];
   real localG[4][4];
   real localXG[4][4];
   real localA[4][4];

   // for non-standart meshes
   //real J[2][2];
   //real Jgrad_i[2];
   //real gradi[2];
   //real Jgrad_j[2];
   //real gradj[2];
   //inline real det_J();
   //real prime_by_var(int what, int varOnFE, int knot_num[4], real ksi, real etta, real tetha);
   //inline int mu(int index);
   //inline int v(int index);
   //inline int nu(int index);
   //real W(int index, real alpha);
   //real d_phi(int index, int what, real ksi, real etta, real tetha);
   //inline real phi(int index, real ksi, real etta, real tetha);
   //void calc_grad(int ij, int index, real ksi, real etta, real tetha);
   //real Integrate2D(const std::function<real(real, real, int, int, int[4])> f, int i, int j, int knot_num[4]);
   //std::function<real(real, real, real, int, int, int[4])> Gij;
   //std::function<real(real, real, real, int, int, int[4])> Mij;

   void AddFirstBounds();
   void AddSecondBounds();
   void AddThirdBounds();
   void AddToGlobalMatricies(element& elem);
   void CreateSLAE(bool isTimed);
   void AssembleMatricies(bool isTimed);
   void CreateM(element& elem);
   void CreateG(element& elem);
   void CreateExtraG(element& elem);
   void Createb(element& elem);

   Matrix* A, *M, *G;
   std::vector<real> b;
   std::vector<real> q1;
   std::vector<real> q2;
   std::vector<real> qt;
   real t_last, th, u0;
   Mesh* mesh;


public:
   int GetKnotsNum() { return num_of_knots; }
   int GetHexasNum() { return num_of_FE; }

   std::vector<real>& GetKnots() { return q1; };
   FEM(Mesh* _mesh);
   void SolveElliptic();
   void SolveParabolic();
   //void GetSolutionOnPlane(real z);
   void Output(std::ofstream& out);

};

