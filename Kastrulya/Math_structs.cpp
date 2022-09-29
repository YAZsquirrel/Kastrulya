#include "Math_structs.h"
#include <set>
#include <array>

namespace maths
{
   Matrix* MakeSparseFormat(int localsize, int elemsize, Mesh* mesh)
   {
      const int N = localsize;
      // set connection table
      std::vector<std::set<int>> map;
      map.resize(elemsize);
      for (auto Elem : mesh->elems)
         for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
               if (Elem.knots_num[i] > Elem.knots_num[j])
                  map[Elem.knots_num[i]].insert(Elem.knots_num[j]);

      Matrix* M = new Matrix;
      M->dim = elemsize;
      M->ig.resize(elemsize + 1, 0);

      for (int i = 0; i < elemsize; i++)
         M->ig[i + 1] = M->ig[i] + map[i].size();
      M->jg.resize(M->ig[elemsize], 0);
      for (int i = 0; i < map.size(); i++)
      {
         std::vector<int> jind;
         jind.reserve(map[i].size());
         std::copy(map[i].begin(), map[i].end(), std::back_inserter(jind));
         for (int j = 0; j < jind.size(); j++)
            M->jg[M->ig[i] + j] = jind[j];
      }

      M->l.resize(M->ig[elemsize], 0.);
      M->u.resize(M->ig[elemsize], 0.);
      M->di.resize(elemsize, 0.);
      return M;

   }

   void AddElement(Matrix* M, int i, int j, real elem)
   {
      bool found = false;
      if (i == j)
         M->di[i] += elem;
      else if (i < j)
      {
         int m;
         for (m = M->ig[j]; m < M->ig[j + 1]; m++)
            if (M->jg[m] == i) { found = true; break; }
         if (found)
            M->u[m] += elem; // i-1?
      }
      else
      {
         int n;
         for (n = M->ig[i]; n < M->ig[i + 1]; n++)
            if (M->jg[n] == j) { found = true; break; }
         if (found)
            M->l[n] += elem; // i-1??
      }
   }

   void MatxVec(std::vector<real>& v, Matrix* M, std::vector<real>& b) // v = M*b
   {
      for (int i = 0; i < M->dim; i++)
         v[i] = M->di[i] * b[i];

      for (int i = 0; i < M->dim; i++)
         for (int j = M->ig[i]; j < M->ig[i + 1]; j++) // -1?
         {
            v[i] += M->l[j] * b[M->jg[j]];
            v[M->jg[j]] += M->u[j] * b[i];
         }
   }

   void copy(std::vector<real>& to, std::vector<real>& from)
   {
      for (int i = 0; i < to.size(); i++)
         to[i] = from[i];
   }
   void copy(std::vector<int>& to, std::vector<int>& from)
   {
      for (int i = 0; i < to.size(); i++)
         to[i] = from[i];
   }

   real scalar(std::vector<real>& v, std::vector<real>& u)
   {
      real sum = 0.;
      for (int i = 0; i < v.size(); i++)
         sum += v[i] * u[i];
      return sum;
   }

   void SolveSLAE(Matrix* M, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> z, r, p, ff, & x = q;
      z.resize(M->dim);
      r.resize(M->dim);
      p.resize(M->dim);
      ff.resize(M->dim);


      real res, alpha, beta, skp, eps = 1e-14;
      int i, k;
      //x = q;

      real lastres;
      MatxVec(ff, M, x);
      for (i = 0; i < M->dim; i++)
         z[i] = r[i] = b[i] - ff[i];
      MatxVec(p, M, z);
      res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));

      for (k = 1; k < 100000 && res > eps; k++)
      {
         lastres = res;
         skp = scalar(p, p);
         alpha = scalar(p, r) / skp;
         for (i = 0; i < M->dim; i++)
         {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         MatxVec(ff, M, r);
         beta = -scalar(p, ff) / skp;
         for (i = 0; i < M->dim; i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = ff[i] + beta * p[i];
         }
         res = sqrt(scalar(r, r)) / sqrt(scalar(b, b));
      }
      //std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void WriteMatrix(Matrix* M)
   {
      double** mat = new double* [M->dim]{};
      for (int i = 0; i < M->dim; i++)
      {
         mat[i] = new double[M->dim]{};
      }

      for (int i = 0; i < M->dim; i++)
      {
         mat[i][i] = M->di[i];
         for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
         {
            mat[i][M->jg[j]] = M->l[j];
            mat[M->jg[j]][i] = M->u[j];
         }
      }

      std::ofstream out("matrix.txt");

      for (int i = 0; i < M->dim; i++)
      {
         for (int j = 0; j < M->dim; j++)
         {
            out.setf(std::ios::left);
            out.width(15);
            out << mat[i][j];
         }
         out << "\n";
      }
   }

   void MatSymmetrisation(Matrix* M, std::vector<real>& b, int i)
   {
      for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
      {
         b[M->jg[j]] -= b[i] * M->u[j];
         M->u[j] = 0;
      }

      for (int n = 0; n < M->dim; n++)
      {
         for (int j = M->ig[n]; j < M->ig[n + 1]; j++)
            if (M->jg[j] == i)
            {
               b[n] -= b[i] * M->l[j];
               M->l[j] = 0.0;
            }
      }
   }
}