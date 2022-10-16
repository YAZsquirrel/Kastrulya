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

   bool operator==(const knot& k) {
      return *this == k;
   }
   bool operator==(const knot* k) {
      bool equal = abs(k->x - x) < 1e-10 && 
                   abs(k->y - y) < 1e-10 &&
                   abs(k->z - z) < 1e-10;

      return equal;
   }
   bool operator!=(const knot& k) {
      return *this != k;
   }

   bool operator!=(const knot* k) {
      bool equal = abs(k->x - x) < 1e-10 &&
         abs(k->y - y) < 1e-10 &&
         abs(k->z - z) < 1e-10;

      return !equal;
   }
};

struct element {
   real lam = 0., gam = 0.;
   const int local_knots_num = 4;
   int n_mat = -1;
   int knots_num[4]{};
   knot* knots[4];

   bool containsKnot(int n)
   {
      bool found = false;
      for (int i = 0; i < local_knots_num && !found; i++)
         found = knots_num[i] == n;
      return found;
   }
   bool containsKnot(knot &k)
   {
      bool found = false;
      for (int i = 0; i < local_knots_num && !found; i++)
         found = knots[i] == &k;
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
   knot& GetVelocity(knot& k, int waterAreaNum)
   {  
      real x1 = areas[waterAreaNum].X1, 
           x2 = areas[waterAreaNum].X2,
           y1 = areas[waterAreaNum].Y1, 
           y2 = areas[waterAreaNum].Y2;

      knot c = knot((x2 + x1) / 2.,
                    (y2 + y1) / 2., 0);
      real klen = sqrt(pow(k.x - c.x, 2.) + pow(k.y - c.y, 2.));
      real len;
      //knot a1 = knot(areas[waterAreaNum].X2 - areas[waterAreaNum].X1,
      //               areas[waterAreaNum].Y2 - areas[waterAreaNum].Y1, 0), 
      //     a2 = a1;
      //     a2.y = -a2.y;

      real yl1 = (k.x - x1) * (x2 - x1) * (y2 - y1) + y1,
           yl2 = (k.x - x1) * (x2 - x1) * (y1 - y2) + y2;

      bool is_lower1 = k.y > yl1,
           is_lower2 = k.y > yl2;
      knot kn = knot(0,0,0);
      real x, y;
      x1 = c.x; x2 = k.x; y1 = c.y; y2 = k.y;
      if (is_lower1)     
         if (is_lower2)
         {
            // v - Верх
            y = areas[waterAreaNum].Y2;
            x = (k.y - y1) * (y - y1) * (x2 - x1) + x1;
            kn.x = -1.;
         }
         else          
         {
            // > - лево
            x = areas[waterAreaNum].X1;
            y = (k.x - x1) * (x - x1) * (y2 - y1) + y1;
            kn.y = -1.;
         }
      else
         if (is_lower2)
         {
            // < - право
            x = areas[waterAreaNum].X2;
            y = (k.x - x1) * (x - x1) * (y2 - y1) + y1;
            kn.y = 1.;
         }
         else
         {
            // ^ - низ
            y = areas[waterAreaNum].Y1;
            x = (k.y - y1) * (y - y1) * (x2 - x1) + x1;
            kn.x = 1.;
         }
      len = sqrt(pow(x - c.x, 2.) + pow(y - c.y, 2.));

      real t = klen / len;
      real v = t > 0.5 ? (1. - t) * 2. * max_v : t * 2. * max_v;
      kn.x *= v;
      kn.y *= v;
      return kn;
   };
   void MakeMesh();

private: 
   void SetBoundConds();
   void SetElemParameters();
   void RemoveNullKnots();
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
   struct boundEdge
   {
      int n_bound;
      knot k1, k2;
      real v1, v2; 
   };

   real max_v;
   std::vector<material> mats;
   std::vector<area> areas;
   std::vector<boundEdge> bes;

   //void Set???();

};

