#pragma once
#include <vector>
typedef double real;

struct bound {
   int knots_num[2];
   real value1;      // ug, th, ub
   real value2;      // beta
};

struct knot
{
   //unsigned int knot_num;
   knot(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
   knot() : x(0.0), y(0.0), z(0.0) {}
   real x, y, z;
   knot& operator=(const knot& k) {
      if (this == &k) return *this;
      x = k.x;
      y = k.y;
      z = k.z;
      return *this;
   }

};

struct element {
   real lam = 0., gam = 0.;
   const int local_knots_num = 4;
   int knots_num[4]{};

   bool containsKnot(int n)
   {
      bool found = false;
      for (int i = 0; i < local_knots_num && !found; i++)
         found = knots_num[i] == n;
      return found;
   }

   element& operator=(const element& elem) {
      if (this == &elem) return *this;
      lam = elem.lam;
      gam = elem.gam;
      for (int i = 0; i < local_knots_num; i++)
         knots_num[i] = elem.knots_num[i];
      return *this;
   }

   element(real gamma, real lambda, int knots_nums[4]) : lam(lambda), gam(gamma)
   {
      for (int i = 0; i < local_knots_num; i++)
         knots_num[i] = knots_nums[i];
   }
   element() : lam(0), gam(0){}
};

class Mesh
{
public:
	std::vector<element> elems;
	std::vector<knot> knots;
	std::vector<bound> bounds1;
	std::vector<bound> bounds2;
	std::vector<bound> bounds3;
   real max_v;
   void MakeMesh();
   void SetBoundConds();
   void SetElemParameters();
   real GetVelocity(knot& k, int waterMatNum)
   {
      
      return...;
   };
private: 
   struct material
   {
      real Cp, Ro, lam, beta;
      material(real _Cp, real _Ro, real _l, real _beta) 
         : Cp(_Cp), Ro(_Ro), lam(_l), beta(_beta){}
   }; 
   struct area
   {
      int n_mat;
      real X1, X2, Y1, Y2;
      area(int nmat, real x1, real x2, real y1, real y2) 
         : n_mat(nmat), X1(x1), X2(x2), Y1(y1), Y2(y2) {}
   };
   std::vector<material> mats;
   std::vector<area> areas;

   //void Set???();

};

