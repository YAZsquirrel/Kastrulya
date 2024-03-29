#pragma once
#include <vector>
typedef double real;

struct bound {
   int knots_num[2];
   int n_mat = -1;
   int n_test = -1;
   real value1;      // ug, th, ub
   real value2;      // beta
};

struct knot
{
   //unsigned int knot_num;
   knot(real _x, real _y, real _z) : x(_x), y(_y), z(_z) {}
   knot(real _x, real _y) : x(_x), y(_y), z(0) {}
   knot() : x(0.0), y(0.0), z(0.0) {}
   real x, y, z;
   knot& operator=(const knot& k) {
      if (this == &k) return *this;
      x = k.x;
      y = k.y;
      z = k.z;
      return *this;
   }

   bool operator==(const knot* k) {
      
      bool equal = (abs(k->x - x) < 1e-10 && 
                   abs(k->y - y) < 1e-10 &&
                   abs(k->z - z) < 1e-10) || *this == k;
      return equal;
   }

   bool operator==(const knot& k) {
      bool equal = (abs(k.x - x) < 1e-10 &&
                    abs(k.y - y) < 1e-10 &&
                    abs(k.z - z) < 1e-10);
      return equal;
   }

   bool operator!=(const knot& k) {
      bool equal = (abs(k.x - x) < 1e-10 &&
                    abs(k.y - y) < 1e-10 &&
                    abs(k.z - z) < 1e-10);
      return !equal;
   }

   bool operator!=(const knot* k) {
      bool equal = (abs(k->x - x) < 1e-10 &&
                    abs(k->y - y) < 1e-10 &&
                    abs(k->z - z) < 1e-10) || *this == k;

      return !equal;
   }
};

struct element {
   real lam = 0., gam = 0.;
   int n_test = 0;
   const int local_knots_num = 4;
   int n_mat = -1;
   int knots_num[4]{};
   knot knots[4];

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
         found = knots[i] == k;
      return found;
   }

   element& operator=(const element& elem) {
      if (this == &elem) return *this;
      lam = elem.lam;
      gam = elem.gam;
      n_test = elem.n_test;
      n_mat = elem.n_mat;
      for (int i = 0; i < local_knots_num; i++)
         knots_num[i] = elem.knots_num[i];
      return *this;
   }

   element(real gamma, real lambda, int _n_mat_test, int knots_nums[4])
      : lam(lambda), gam(gamma), n_test(_n_mat_test)
   {
      for (int i = 0; i < local_knots_num; i++)
         knots_num[i] = knots_nums[i];
   }
   element() : lam(0), gam(0), n_test(0){}
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

      knot kn = knot(0.0, 0.0, 0.0);
      if (x2 <= k.x || k.x <= x1 || y2 <= k.y || k.y <= y1)
         return kn;

      knot c = knot((x2 + x1) / 2.,
                    (y2 + y1) / 2., 0);
      real klen = sqrt(pow(k.x - c.x, 2.) + pow(k.y - c.y, 2.));
      real len;
      //knot a1 = knot(areas[waterAreaNum].X2 - areas[waterAreaNum].X1,
      //               areas[waterAreaNum].Y2 - areas[waterAreaNum].Y1, 0), 
      //     a2 = a1;
      //     a2.y = -a2.y;

      real yl1 = (k.x - x1) / (x2 - x1) * (y2 - y1) + y1,
           yl2 = (k.x - x1) / (x2 - x1) * (y1 - y2) + y2;

      bool is_lower1 = k.y > yl1,
           is_lower2 = k.y > yl2;
      real x, y;
      x1 = c.x; x2 = k.x; y1 = c.y; y2 = k.y;
      if (is_lower1)     
         if (is_lower2)
         {
            // v - ����
            y = areas[waterAreaNum].Y2;
            //x = (k.y - y1) / (y - y1) * (x2 - x1) + x1;
            x = (y - y1) / (y2 - y1) * (x2 - x1) + x1;
            kn.x = -1.;
         }
         else          
         {
            // > - ����
            x = areas[waterAreaNum].X1;
            //y = (k.x - x1) / (x - x1) * (y2 - y1) + y1;
            y = (x - x1) / (x2 - x1) * (y2 - y1) + y1;
            kn.y = -1.;
         }
      else
         if (is_lower2)
         {
            // < - �����
            x = areas[waterAreaNum].X2;
            //y = (k.x - x1) / (x - x1) * (y2 - y1) + y1;
            y = (x - x1) / (x2 - x1) * (y2 - y1) + y1;
            kn.y = 1.;
         }
         else
         {
            // ^ - ���
            y = areas[waterAreaNum].Y1;
            //x = (k.y - y1) / (y - y1) * (x2 - x1) + x1;
            x = (y - y1) / (y2 - y1) * (x2 - x1) + x1;
            kn.x = 1.;
         }
      len = sqrt(pow(x - c.x, 2.) + pow(y - c.y, 2.));

      real t = klen / len;
      real v = t > 0.5 ? (1. - t) * 2. * max_v : t * 2. * max_v;
      //real v = (1. - t) * max_v;
      kn.x *= v;
      kn.y *= v;
      return kn;
   };
   bool isInTriangle(knot& k1, knot& k2, knot& k3, knot& p)
   {
      knot k12 = knot(k2.x - k1.x, k2.y - k1.y),
         k23 = knot(k3.x - k2.x, k3.y - k2.y),
         k31 = knot(k1.x - k3.x, k1.y - k3.y),
         pk1 = knot(p.x - k1.x, p.y - k1.y),
         pk2 = knot(p.x - k2.x, p.y - k2.y),
         pk3 = knot(p.x - k3.x, p.y - k3.y);
      real pk12 = Cross(pk1, k12).z, pk23 = Cross(pk2, k23).z, pk31 = Cross(pk3, k31).z;

      return (pk12 < 1e-12 && pk23 < 1e-12 && pk31 < 1e-12) || (pk12 >= 1e-12 && pk23 >= 1e-12 && pk31 >= 1e-12);
   }
   knot& Cross(knot& k1, knot& k2)
   {
      knot k = knot(k1.y * k2.z - k1.z * k2.y, k1.x * k2.z - k1.z * k2.x, k1.x * k2.y - k1.y * k2.x);
      return k;
   }
   knot& intersection(knot p1, knot p2, knot q1, knot q2)
   {
      knot k;
      real d = (p1.x - p2.x) * (q1.y - q2.y) - (q1.x - q2.x) * (p1.y - p2.y);
      k.x = ((p1.x * p2.y - p1.y * p2.x) * (q1.x - q2.x) - (q1.x * q2.y - q1.y * q2.x) * (p1.x - p2.x)) / d;
      k.y = ((p1.x * p2.y - p1.y * p2.x) * (q1.y - q2.y) - (q1.x * q2.y - q1.y * q2.x) * (p1.y - p2.y)) / d;

      return k;
   }

   void MakeMesh();

private: 
   void SetBoundConds();
   void SetElemParameters();
   void RemoveNullKnots();
   struct material
   {
      real Cp, Ro, lam;
      int n_test;
      material(real _Cp, real _Ro, real _l,  int _n_test) 
         : Cp(_Cp), Ro(_Ro), lam(_l), n_test(_n_test) {}
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

   real max_v = 0;
   std::vector<material> mats;
   std::vector<area> areas;
   std::vector<boundEdge> bes;
   std::vector<real> X, Y;

   //void Set???();

};

