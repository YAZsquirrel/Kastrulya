#include "Math_structs.h"
#include <set>
#include <array>
#include <cmath>
#include <ppl.h>
#include <time.h>

using namespace concurrency;


namespace maths
{
   Matrix* MakeSparseFormat(int localsize, int size, Mesh* mesh)
   {
      const int N = localsize;
      // set connection table
      std::vector<std::set<int>> map;
      map.resize(size);
      for (auto& Elem : mesh->elems)
         for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
               if (Elem.knots_num[i] > Elem.knots_num[j])
                  map[Elem.knots_num[i]].insert(Elem.knots_num[j]);

      Matrix* M = new Matrix;
      M->dim = size;
      M->ig.resize(size + 1, 0);

      for (int i = 0; i < size; i++)
         M->ig[i + 1] = M->ig[i] + map[i].size();
      M->jg.resize(M->ig[size], 0);
      for (int i = 0; i < map.size(); i++)
      {
         std::vector<int> jind;
         jind.reserve(map[i].size());
         std::copy(map[i].begin(), map[i].end(), std::back_inserter(jind));
         for (int j = 0; j < jind.size(); j++)
            M->jg[M->ig[i] + j] = jind[j];
      }

      M->l.resize(M->ig[size], 0.);
      M->u.resize(M->ig[size], 0.);
      M->di.resize(size, 0.);
      M->isDense = false;
      return M;

   }

   Matrix* MakeDenseFormat(int size)
   {
      Matrix* A = new Matrix();
      A->dim = size;
      A->dense.resize(size);
      for (size_t i = 0; i < A->dim; i++)
         A->dense[i].resize(A->dim, 0.);
      return A;
   }

   void AddElement(Matrix* M, int i, int j, real elem)
   {
      if (!M->isDense)
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
      else
         M->dense[i][j] += elem;
   }

   void MatxVec(std::vector<real>& v, Matrix* M, std::vector<real>& b) // v = M*b
   {
      if (!M->isDense)
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
      else
      {
         for (size_t i = 0; i < M->dim; i++)
         {
            real sum = 0;
            #pragma omp parallel for reduction(+:sum)
            for (int j = 0; j < M->dim; j++)
               sum += M->dense[i][j] * b[j];
            v[i] = sum;
         }
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

      real lastres;
      MatxVec(Ar, M, x);
      for (i = 0; i < M->dim; i++)
         z[i] = r[i] = b[i] - Ar[i];
      MatxVec(p, M, z);
      real b_norm = sqrt(scalar(b, b));
      res = sqrt(scalar(r, r)) / b_norm;

      for (k = 1; k < 100000 && res > eps; k++)
      {
         skp = scalar(p, p);
         alpha = scalar(p, r) / skp;

#pragma omp parallel for
         for (i = 0; i < M->dim; i++)
         {
            x[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         MatxVec(Ar, M, r);
         beta = -scalar(p, Ar) / skp;
#pragma omp parallel for
         for (i = 0; i < M->dim; i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = Ar[i] + beta * p[i];
         }
         res = sqrt(scalar(r, r)) / b_norm;
      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void SolveSLAE_LU(Matrix* A, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> y, &x = q;
      y.resize(A->dim, 0);
      std::clock_t start = clock();

      Matrix *LU = MakeDenseFormat(A->dim);
      // L 

//#pragma omp parallel for
      for (int i = 0; i < LU->dim; i++)
      {
//#pragma omp parallel for
         for (int j = 0; j < LU->dim; j++)
         {
            if (i <= j) // U
            {
               real sum = 0.;
#pragma omp parallel for reduction (+:sum)
               for (int k = 0; k < i; k++)
                  sum += LU->dense[i][k] * LU->dense[k][j];
               LU->dense[i][j] = A->dense[i][j] - sum;
            }
            else // L
            {
               real sum = 0.;
#pragma omp parallel for reduction (+:sum)
               for (int k = 0; k < j; k++)
                  sum += LU->dense[i][k] * LU->dense[k][j];
               LU->dense[i][j] = (A->dense[i][j] - sum) / LU->dense[j][j];
            }
         }
      }
      SolveForL(y, b, LU);
      SolveForU(x, y, LU);

      std::clock_t end = clock();

      real res = 0;
      y.resize(A->dim, 0);
      MatxVec(y, A, x);
      for (size_t i = 0; i < A->dim; i++)
         y[i] -= b[i];
      res = sqrt(scalar(y, y) / scalar(b, b));
      std::cout << res << '\n';
      std::cout << "Work time (N = " << A->dim << "): " << end - start << "ms [" << (end - start) / 1000 << "s] (" << (end - start) / 60000 << "m)\n";

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

      SolveForL(r, t2, SQ);
      SolveForU(z, r, SQ);

      MatxVec(t1, A, z);                            
      SolveForL(p, t1, SQ);

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
         SolveForL(t1, r, SQ);                        // t2 = Q1*rk -> Q*t2 = rk
         MatxVec(At1, A, t1);                       // A*t2 = At
         SolveForU(t2, At1, SQ);                      // t1 = S1*At -> S*t1 = At

         beta = -scalar(t2, p) / skp;              // beta = (pk0, t) / (pk0, pk0)          

         SolveForU(t2, r, SQ);                        // Q*zk1 = rk1
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
      if (!M->isDense)
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
      else 
      {
         for (int j = 0; j < M->dim; j++)
         {
            b[j] -= b[i] * M->dense[i][j];
            M->dense[i][j] = 0;
         }
         M->dense[i][i] = 1.;
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
      M->isDense = false;

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

   void SolveForL(std::vector<real>& q, std::vector<real>& b, Matrix* M) // y = 1/L * b
   {
      if (!M->isDense)
      {
         for (int k = 1, k1 = 0; k < M->dim; k++, k1++)
         {
            double sum = 0;
            for (int i = M->ig[k1]; i < M->ig[k]; i++)
               sum += M->l[i] * q[M->jg[i]];
   
            q[k1] = b[k1] - sum;
         }
      }
      else
      {
//#pragma omp parallel for  
         for (int i = 0; i < M->dim; i++)
         {  
            real sum = 0;
#pragma omp parallel for reduction (+:sum)
            for (int j = 0; j < i; j++)
               sum += q[j] * M->dense[i][j];

            q[i] = b[i] - sum;
         }
      }
   }
   void SolveForU(std::vector<real>& q, std::vector<real>& b, Matrix* M) // x = 1/U * y
   {
      if (!M->isDense)
      {
         for (int k = M->dim - 1; k > 0; k--)
         {
            /// ??????????????????????????????????????
            //real sum = b[k1];
            ////q[k1] = b[k1] / M->di[k1];
            ////double v_el = q[k1];
            //
            //for (int i = M->ig[k1]; i < M->ig[k]; i++)
            //   b[M->jg[i]] -= M->l[i] * v_el;

            real sum = b[k];

            for (int j = M->ig[k]; j < M->ig[k + 1]; j++) // -1?
               sum -= M->l[j] * q[M->jg[j]];
            
            q[k] = sum / M->di[k];
         } 
      }
      else
      {
//#pragma omp parallel for reduction(-:sum)
         for (int i = M->dim - 1; i >= 0; i--)
         {
            real sum = b[i];
#pragma omp parallel for reduction(-:sum)
            for (int j = i + 1; j < M->dim; j++)
               sum -= q[j] * M->dense[i][j];

            q[i] = sum / M->dense[i][i];
         }
      }
   }
}