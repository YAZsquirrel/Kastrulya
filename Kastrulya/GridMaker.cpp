#include "GridMaker.h"
#include <fstream>

void Mesh::MakeMesh()
{
   std::ifstream fmat("Materials.txt");
   int k;
   fmat >> k;
   mats.reserve(k);

   for (int i = 0; i < k; i++)
   {
      real Cp, Ro, l, v;
      int n_test;
      fmat >> Cp >> Ro >> l  >> n_test >> v;
      if (v > 1e-4)
         max_v = v;

      material mat(Cp, Ro, l, n_test);
      mats.push_back(mat);
   }
   fmat.close();

   std::ifstream fgrid("GridDescription.txt");
   int nX, nY;
   fgrid >> nX >> nY;

   std::vector<real> Xs, Ys;
   std::vector<real> Xr, Yr;
   X.reserve(nX);
   Y.reserve(nY);
   Xs.reserve(nX - 1);
   Ys.reserve(nY - 1);
   Xr.reserve(nX - 1);
   Yr.reserve(nY - 1);

   for (size_t i = 0; i < nX; i++)
   {
      real x;
      fgrid >> x;
      X.push_back(x);
   }
   for (size_t i = 0; i < nY; i++)
   {
      real y;
      fgrid >> y;
      Y.push_back(y);
   }
   for (size_t i = 0; i < nX - 1; i++)
   {
      real xs;
      fgrid >> xs;
      Xs.push_back(xs);
   }
   for (size_t i = 0; i < nY - 1; i++)
   {
      real ys;
      fgrid >> ys;
      Ys.push_back(ys);
   }
   for (size_t i = 0; i < nX - 1; i++)
   {
      real xr;
      fgrid >> xr;
      Xr.push_back(xr);
   }
   for (size_t i = 0; i < nY - 1; i++)
   {
      real yr;
      fgrid >> yr;
      Yr.push_back(yr);
   }

   int nBound;
   fgrid >> nBound;
   for (size_t i = 0; i < nBound; i++)
   {
      int bound_n;
      real value1, value2 = 0.;
      int p1, p2, q;
      bool axis;
      fgrid >> bound_n;
      if (bound_n == 3)
         fgrid >> value1 >> value2 >> p1 >> p2 >> q >> axis;
      else
         fgrid >> value1 >> p1 >> p2 >> q >> axis;
         
      boundEdge b;
      if (!axis)
      {
         b.k1 = knot(X[p1], Y[q], 0.0);
         b.k2 = knot(X[p2], Y[q], 0.0);
      }
      else
      {
         b.k1 = knot(X[q], Y[p1], 0.0);
         b.k2 = knot(X[q], Y[p2], 0.0);
      }
      b.v1 = value1;
      b.v2 = value2;
      b.n_bound = bound_n;
      bes.push_back(b);
   }


   int nAreas;
   fgrid >> nAreas;

   for (size_t i = 0; i < nAreas; i++)
   {
      int mat_n; 
      int x1, x2, y1, y2;
      fgrid >> mat_n >> x1 >> x2 >> y1 >> y2;
      area a(mat_n, X[x1], X[x2], Y[y1], Y[y2]);
      areas.push_back(a);
   }
   fgrid.close();

   //код для создания сетки
   knots.reserve(nX * nY);
   elems.reserve((nX - 1) * (nY - 1));

   int x_size = X.size(), y_size = Y.size();
   for (size_t i = 0; i < Xs.size(); i++)
      x_size += Xs[i] - 1;
   for (size_t i = 0; i < Ys.size(); i++)
      y_size += Ys[i] - 1;
   
   std::vector<real> X_full, Y_full;
   X_full.reserve(x_size);
   Y_full.reserve(y_size);

   X_full.push_back(X[0]);
   for (size_t i = 1; i < X.size(); i++)
   {
      bool isRegular = !(abs(Xr[i - 1] - 1.) > 1e-10);
      real xl = X[i - 1],
           xr = X[i],
           xh = (xr - xl) / Xs[i - 1];
      real xb1 = isRegular ? xh :
                             xh * (1. - Xr[i - 1]) / (1. - pow(Xr[i - 1], Xs[i - 1]));
      real x = xl + xb1;
      X_full.push_back(x);
      for (size_t n = 1; n < Xs[i - 1]; n++)
      {
         real scale = isRegular ? 1. : pow(Xr[i - 1], n);
         x += xb1 * scale;
         X_full.push_back(x);
      }
   }

   Y_full.push_back(Y[0]);
   for (size_t j = 1; j < Y.size(); j++)
   {
      bool isRegular = !(abs(Yr[j - 1] - 1.) > 1e-10);
      real yl = Y[j - 1], 
           yr = Y[j],
           yh = (yr - yl) / Ys[j - 1];
      real yb1 = isRegular ? yh :
                             yh * (1. - Yr[j - 1]) / (1. - pow(Yr[j - 1], Ys[j - 1]));
      real y = yl + yb1;
      Y_full.push_back(y);
      for (size_t n = 1; n < Ys[j - 1]; n++)
      {
         real scale = isRegular ? 1. : pow(Yr[j - 1], n);
         y += yb1 * scale;
         Y_full.push_back(y);
      }
   }

   for (size_t iy = 0; iy < y_size; iy++)
      for (size_t ix = 0; ix < x_size; ix++)
      {
         knot kn = knot(X_full[ix], Y_full[iy], 0);
         knots.push_back(kn);
      }

   //std::ofstream of("knots.txt");
   //for (size_t i = 0; i < knots.size(); i++)
   //{
   //   of << knots[i].x << " " << knots[i].y << '\n';
   //}
   //of.close();

   for (size_t j = 0; j < y_size - 1; j++)
      for (size_t i = 0; i < x_size - 1; i++)
      {
         int n_mat = -1;
         for (auto& a : areas)
         {
            bool found = false;
            int ks[4]{i + j * x_size, (i + 1) + j * x_size,  i + (j + 1) * x_size, (i + 1) + (j + 1) * x_size };

            real xm = (knots[ks[0]].x + knots[ks[1]].x) / 2., 
                 ym = (knots[ks[0]].y + knots[ks[2]].y) / 2.;

            for (size_t i = 0; i < 4 && !found; i++)
               found = a.X1 < xm && xm < a.X2 &&
                       a.Y1 < ym && ym < a.Y2;
            
            if (found) 
            { 
               n_mat = a.n_mat; 
               break;
            }
         }
         if (n_mat < 0)
            continue;
         element e;
         e.knots_num[0] = i       + j       * x_size;
         e.knots_num[1] = (i + 1) + j       * x_size;
         e.knots_num[2] = i       + (j + 1) * x_size;
         e.knots_num[3] = (i + 1) + (j + 1) * x_size;
         elems.push_back(e);
      }

   for (auto& el : elems)
      for (size_t i = 0; i < 4; i++)
         el.knots[i] = knots[el.knots_num[i]];

   SetElemParameters();
   RemoveNullKnots();
   SetBoundConds();
}

void Mesh::SetBoundConds()
{
   for (auto& e : elems)
   {
      for (auto& be : bes)
      {
         knot ks[4] = {knots[e.knots_num[0]], knots[e.knots_num[1]], 
                       knots[e.knots_num[2]], knots[e.knots_num[3]] };
         int kn[2]{};

         bool onEdge = false;
         if (abs(be.k1.x - be.k2.x) > 1e-10) // ------
         {  
            if (onEdge = abs(be.k1.y - ks[0].y) < 1e-10 // на том же yl
               && (be.k1.x - ks[0].x < 1e-10 && ks[0].x - be.k2.x < 1e-10)
               && (be.k1.x - ks[1].x < 1e-10 && ks[1].x - be.k2.x < 1e-10))
            {
               kn[0] = e.knots_num[0];
               kn[1] = e.knots_num[1];
            }
            else if (onEdge = abs(be.k1.y - ks[2].y) < 1e-10 // на том же yr
               && (be.k1.x - ks[2].x < 1e-10 && ks[2].x - be.k2.x < 1e-10)
               && (be.k1.x - ks[3].x < 1e-10 && ks[3].x - be.k2.x < 1e-10))
            {
               kn[0] = e.knots_num[2];
               kn[1] = e.knots_num[3];
            }
         }
         else
         {
            if (onEdge = abs(be.k1.x - ks[0].x) < 1e-10 // на том же xl
               && (be.k1.y - ks[0].y < 1e-10 && ks[0].y - be.k2.y < 1e-10)
               && (be.k1.y - ks[2].y < 1e-10 && ks[2].y - be.k2.y < 1e-10))
            {
               kn[0] = e.knots_num[0];
               kn[1] = e.knots_num[2];
            }
            else if (onEdge = abs(be.k1.x - ks[1].x) < 1e-10 // на том же xr
               && (be.k1.y - ks[1].y < 1e-10 && ks[1].y - be.k2.y < 1e-10)
               && (be.k1.y - ks[3].y < 1e-10 && ks[3].y - be.k2.y < 1e-10))
            {
               kn[0] = e.knots_num[1];
               kn[1] = e.knots_num[3];
            }
         }
         if (onEdge)
         {
            bound b;
            b.knots_num[0] = kn[0];
            b.knots_num[1] = kn[1];
            b.value1 = be.v1;
            b.value2 = be.v2;
            b.n_mat = e.n_mat;
            b.n_test = e.n_test;
            switch (be.n_bound)
            {
               case 1: 
                  bounds1.push_back(b); 
                  break;
               case 2: 
                  bounds2.push_back(b); 
                  break;
               case 3: 
                  bounds3.push_back(b);
                  break;
            }
         }
      }
   }

}

void Mesh::SetElemParameters()
{
   for (auto& e : elems)
   {
      int n_mat = -1;
      knot center = knot((e.knots[0].x + e.knots[3].x) / 2., 
                              (e.knots[0].y + e.knots[3].y) / 2., 0);
      for (auto& a : areas)
      {
         if (a.X1 < center.x && center.x < a.X2 &&
             a.Y1 < center.y && center.y < a.Y2)
             {n_mat = a.n_mat; break;}
         else n_mat = -1;
      }
      if (n_mat > -1)
      {
         e.n_mat = n_mat;
         e.gam = mats[n_mat].Cp * mats[n_mat].Ro;
         e.lam = mats[n_mat].lam;
         e.n_test = mats[n_mat].n_test;
      }
   }

}

void Mesh::RemoveNullKnots()
{
   //if there's knots that are orphans
   std::vector<int> toRemove;
   toRemove.reserve(knots.size() / 2);

   for (auto& k : knots)
   {
      bool isOrphan = true;
      for (auto& el : elems)
         if (el.containsKnot(k))
         {
            isOrphan = false;
            break;
         }

      for (int i = 0; i < knots.size() && isOrphan; i++)
         if (&knots[i] == &k)
         {
            toRemove.push_back(i);
            break;
         }
   }

   //remove 'em all
   for (size_t j = 0; j < toRemove.size(); j++)
   {
      knots.erase(knots.begin() + toRemove[j]);
      for (size_t jr = j; jr < toRemove.size(); jr++)
         toRemove[jr]--;

   }

   //renumerate
   for (auto& el : elems)
      for (size_t k = 0; k < 4; k++)
         if (el.knots_num[k] >= knots.size())
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == el.knots[k])
               {
                  el.knots_num[k] = i;
                     break;
               }
         }
         else if (knots[el.knots_num[k]] != el.knots[k])
         {
            for (size_t i = 0; i < knots.size(); i++)
               if (knots[i] == el.knots[k])
               {
                  el.knots_num[k] = i;
                  break;
               }    
         }
}
