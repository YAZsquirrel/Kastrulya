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
      real Cp, Ro, l, beta, v;
      fmat >> Cp >> Ro >> l >> beta >> v;
      if (v > 1e-10)
         max_v = v;

      material mat(Cp, Ro, l, beta);
      mats.push_back(mat);
   }
   fmat.close();

   std::ifstream fgrid("GridDescription.txt");
   int nX, nY;
   fgrid >> nX >> nY;

   std::vector<real> X, Y;
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

   int nAreas;
   fgrid >> nAreas;

   for (size_t i = 0; i < nAreas; i++)
   {
      int mat_n; 
      real x1, x2, y1, y2;
      fgrid >> mat_n >> x1 >> x2 >> y1 >> y2;
      area a(mat_n, x1, x2, y1, y2);
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
   X_full.resize(x_size);
   Y_full.resize(y_size);

   X_full.push_back(X[0]);
   for (size_t i = 1; i < X.size(); i++)
   {
      real xl = X[i - 1],
           xr = X[i],
           xh = (xr - xl) / Xs[i - 1];
      real xb1 = xh * (1. - Xr[i - 1]) / (1. - pow(Xr[i - 1], Xs[i - 1]));
      real x = xl + xb1;
      for (size_t n = 1; n < Xs[i - 1]; n++)
      {
         x += xb1 * pow(Xr[i - 1], n);
         X_full.push_back(x);
      }
   }

   Y_full.push_back(Y[0]);
   for (size_t j = 1; j < Y.size(); j++)
   {
      real yl = Y[j - 1], 
           yr = Y[j],
           yh = (yr - yl) / Ys[j - 1];
      real yb1 = yh * (1. - Yr[j - 1]) / (1. - pow(Yr[j - 1], Ys[j - 1]));
      real y = yl + yb1;
      for (size_t n = 1; n < Ys[j - 1]; n++)
      {
         y += yb1 * pow(Yr[j - 1], n);
         Y_full.push_back(y);
      }
      Y_full.push_back(Y[j]);
   }

   for (size_t iy = 0; iy < y_size; iy++)
      for (size_t ix = 0; ix < x_size; ix++)
      {
         knot kn(X_full[ix], Y_full[iy], 0);
         knots.push_back(kn);
      }

   for (size_t j = 0; j < nY - 1; j++)
      for (size_t i = 0; i < nX - 1; i++)
      {
         element e;
         e.knots_num[0] = i       + j       * nX;
         e.knots_num[1] = (i + 1) + j       * nX;
         e.knots_num[2] = i       + (j + 1) * nX;
         e.knots_num[3] = (i + 1) + (j + 1) * nX;
         elems.push_back(e);
      }



}

void Mesh::SetBoundConds()
{
}

void Mesh::SetElemParameters()
{
   for (size_t e = 0; e < elems.size(); e++)
   {
      
   }

}
