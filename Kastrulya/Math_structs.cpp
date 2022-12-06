#include "Math_structs.h"
#include <set>
#include <array>
#include <cmath>

namespace maths
{
   Matrix* MakeSparseFormat(int localsize, int elemsize, Mesh* mesh)
   {
      const int N = localsize;
      // set connection table
      std::vector<std::set<int>> map;
      map.resize(elemsize);
      for (auto& Elem : mesh->elems)
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

   void SolveSLAE_LOS(Matrix* M, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> z, r, p, Ar, &x = q;
      z.resize(M->dim);
      r.resize(M->dim);
      p.resize(M->dim);
      Ar.resize(M->dim);


      real res, alpha, beta, skp, eps = 1e-14;
      int i, k;
      //x = q;

      SolveSLAE_Relax(M, x, b, 1.6);
      real lastres;
      MatxVec(Ar, M, x);
      for (i = 0; i < M->dim; i++)
         z[i] = r[i] = b[i] - Ar[i];
      MatxVec(p, M, z);
      real b_norm = sqrt(scalar(b, b));
      res = sqrt(scalar(r, r)) / b_norm;

      for (k = 1; k < 1000000 && res > eps; k++)
      {
         skp = scalar(p, p);
         alpha = scalar(p, r) / skp;
         for (i = 0; i < M->dim; i++)
         {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         MatxVec(Ar, M, r);
         beta = -scalar(p, Ar) / skp;
         for (i = 0; i < M->dim; i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = Ar[i] + beta * p[i];
         }
         res = sqrt(scalar(r, r)) / b_norm;
         //SolveSLAErelax(M, x, b, 1.6);
      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void SolveSLAE_LOSnKholessky(Matrix* A, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> z, r, p, At1, t2, t1, &x = q;
      z.resize(A->dim);
      r.resize(A->dim);
      p.resize(A->dim);
      At1.resize(A->dim);
      t2.resize(A->dim);
      t1.resize(A->dim);
      real res, alpha, beta, skp, eps = 1e-14;
      int i, k;
      //x = q;

      Matrix* SQ = MakeKholessky(A);

      MatxVec(t1, A, x);                            
      real b_norm;

      for (int i = 0; i < A->dim; i++)         
         t2[i] = b[i] - t1[i];

      SolveL(r, t2, SQ);
      SolveLT(z, r, SQ);

      MatxVec(t1, A, z);                            
      SolveL(p, t1, SQ);

      b_norm = sqrt(scalar(b, b));
      res = sqrt(scalar(r, r)) / b_norm;

      for (k = 1; k < 100000 && res > eps; k++)
      {
         skp = scalar(p, p);
         alpha = scalar(p, r) / skp;               // a = (pk0, rk0) / (pk0, pk0)

         for (int i = 0; i < A->dim; i++)          // xk1 = xk0 + a*zk0
         {                                         // rk1 = rk0 - a*pk0 
            x[i] += alpha * z[i];                     
            r[i] -= alpha * p[i];                    
         }

         //MatxVec(Ar, M, r);                        
         //SolveLLT(t1, Ar, SQ);                   //// S1*A*Q1*rk = t1 =>        
         SolveL(t1, r, SQ);                        // t2 = Q1*rk -> Q*t2 = rk
         MatxVec(At1, A, t1);                       // A*t2 = At
         SolveLT(t2, At1, SQ);                      // t1 = S1*At -> S*t1 = At

         beta = -scalar(t2, p) / skp;              // beta = (pk0, t) / (pk0, pk0)          

         SolveLT(t2, r, SQ);                        // Q*zk1 = rk1
         for (int i = 0; i < A->dim; i++)
         {
            z[i] = t1[i] + beta * z[i];                   // zk1 += beta
            p[i] = t2[i] + beta * p[i];
         }
         
         res = sqrt(scalar(r, r)) / b_norm;
         if (res != res)
            std::cerr << "Error: NaN detected!" << std::endl;

      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
      delete SQ;
   }


   void SolveSLAE_Relax(Matrix* M, std::vector<real>& q, std::vector<real>& b, real w)
   {
      w = w > 2. ? 2. : w; 
      w = w <= 0. ? 0.5 : w; 
      std::vector<real> Ar, q1, q2;
      Ar.resize(M->dim);
      q1.resize(M->dim);
      q2.resize(M->dim);

      MatxVec(Ar, M, q);
      copy(q1, q);

      real res = 0., prevres = 1e7;
      for (size_t i = 0; i < M->dim; i++)
         res += b[i] - Ar[i];
      int k = 0;

      for (; res - prevres < 1e-14 && k < 100000; k++)
      {         
         copy(q2, q1);
         for (size_t i = 0; i < M->dim; i++)
         {
            Ar[i] = M->di[i] * q1[i];
            for (int j = M->ig[i]; j < M->ig[i + 1]; j++) // -1?
            {
               Ar[i] += M->l[j] * q1[M->jg[j]];
               Ar[M->jg[j]] += M->u[j] * q1[i];
            }

            q1[i] += w / M->di[i] * (b[i] - Ar[i]);
         }
         prevres = res;
         res = 0.;
         MatxVec(Ar, M, q1);

         for (size_t i = 0; i < M->dim; i++)
            res += b[i] - Ar[i];
      }
      copy(q, q2);
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

   Matrix* MakeKholessky(Matrix *A) 
   {
      double sum_d, sum_l;
      Matrix* M = new Matrix();
      M->dim = A->dim;
      M->ig.resize(A->ig.size());
      M->jg.resize(A->jg.size());
      M->l.resize(A->l.size());
      M->u.resize(A->u.size());
      M->di.resize(A->di.size());

      copy(M->di, A->di);
      copy(M->ig, A->ig);
      copy(M->jg, A->jg);
      copy(M->u, A->u);
      copy(M->l, A->l);
      

      for (int k = 0; k < A->dim; k++)
      {
         sum_d = 0;
         int i_s = M->ig[k], i_e = M->ig[k + 1];

         for (int i = i_s; i < i_e; i++)
         {
            sum_l = 0;
            int j_s = M->ig[M->jg[i]], j_e = M->ig[M->jg[i] + 1];
            for (int m = i_s; m < i; m++)
               for (int j = j_s; j < j_e; j++)
                  if (M->jg[m] == M->jg[j])
                  {
                     sum_l += M->l[m] * M->l[j];
                     j_s++;
                  }
            M->l[i] = (M->l[i] - sum_l) / M->di[M->jg[i]];
            sum_d += M->l[i] * M->l[i];
         }
         M->di[k] = sqrt(abs(M->di[k] - sum_d));

      }
      return M;
   }
   void SolveLLT(std::vector<real>& q, std::vector<real>& b, Matrix* M)
   {
      SolveL(b, q, M);
      SolveLT(q, q, M);
   }
   void SolveL(std::vector<real>& q, std::vector<real>& b, Matrix* M)
   {
      for (int k = 1, k1 = 0; k <= M->dim; k++, k1++)
      {
         double sum = 0;

         for (int i = M->ig[k1]; i < M->ig[k]; i++)
            sum += M->l[i] * q[M->jg[i]];

         q[k1] = (b[k1] - sum) / M->di[k1];
      }
   }
   void SolveLT(std::vector<real>& q, std::vector<real>& b, Matrix* M)
   {
      for (int k = M->dim, k1 = M->dim - 1; k > 0; k--, k1--)
      {

         q[k1] = b[k1] / M->di[k1];
         double v_el = q[k1];

         for (int i = M->ig[k1]; i < M->ig[k]; i++)
            b[M->jg[i]] -= M->l[i] * v_el;
      }
   }
}