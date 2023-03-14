#include "Math_structs.h"
#include <set>
#include <array>
#include <cmath>
#include <atomic>
#include <numeric>
#include <array>
#include <time.h>


namespace maths
{
   Matrix* MakeSparseRowColumnFormat(int localsize, int size, Mesh* mesh)
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
      M->format = MatrixFormat::SparseRowColumn;
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

   void MakeSparseRowFormatFromRCF(Matrix* M, Matrix*& M_srf)
   {
       Matrix* SRF = M_srf = new Matrix();

       SRF->dim = M->dim;
       SRF->format = MatrixFormat::SparseRow;

       SRF->ig.resize(SRF->dim + 1);
       SRF->jg.resize(M->ig[SRF->dim] + SRF->dim);
       SRF->gg.resize(M->ig[SRF->dim] + SRF->dim, 0);
       copy(SRF->di, M->di);

       std::ofstream logfile;
       logfile.open("pardiso64.log");
       if (!logfile)
       {
          std::cerr << "Cannot open pardiso64.log\n";
          return;
       }

       ConvertFromSRCFToSRF(M, SRF);
   }

   void ConvertFromSRCFToSRF(Matrix* M_srcf, Matrix* M_srf) // rsf -> g, csr -> a
   {
       std::vector<int> adr;
       int N = M_srcf->dim;
       // подсчитываем число элементов в каждой строчке
       adr.resize(N, 0);

       for (int i = 0; i < N; i++)
       {
           adr[i] += 1; // диагональ
           // верхний треугольник
           for (int j = M_srcf->ig[i]; j < M_srcf->ig[i + 1]; j++)
               adr[M_srcf->jg[j]]++;
       }
       // ia
       M_srf->ig[0] = 0;
       for (int i = 0; i < N; i++)
          M_srf->ig[i + 1] = M_srf->ig[i] + adr[i];

       // ja,  a
       for (int i = 0; i < M_srcf->ig[N] + N; i++)
          M_srf->gg[i] = 0;

       for (int i = 0; i < N; i++)
           adr[i] = M_srf->ig[i]; // в какую позицию заносить значение

       // диагональ
       for (int i = 0; i < N; i++)
       {
          M_srf->jg[adr[i]] = i;
          M_srf->gg[adr[i]] = M_srcf->di[i];
           adr[i]++;
       }

       // верхний треугольник
       for (int i = 0; i < N; i++)
       {
           for (int j = M_srcf->ig[i]; j < M_srcf->ig[i + 1]; j++)
           {
               int k = M_srcf->jg[j];
               M_srf->jg[adr[k]] = i;
               M_srf->gg[adr[k]] = M_srcf->u[j];
               
               adr[k]++;
           }
       }
   }

   void EqualizeRSFToCSR(Matrix* M_srcf, Matrix* M_srf)
   {
      std::vector<int> adr;
      int N = M_srcf->dim;
      // подсчитываем число элементов в каждой строчке
      adr.resize(N, 0);

      for (int i = 0; i < N; i++)
      {
         adr[i] += 1; // диагональ
         // верхний треугольник
         for (int j = M_srcf->ig[i]; j < M_srcf->ig[i + 1]; j++)
            adr[M_srcf->jg[j]]++;
      }

      for (int i = 0; i < N; i++)
         adr[i] = M_srf->ig[i]; // в какую позицию заносить значение

      // диагональ
      for (int i = 0; i < N; i++)
      {
         M_srf->gg[adr[i]] = M_srcf->di[i];
         adr[i]++;
      }

      // верхний треугольник
      for (int i = 0; i < N; i++)
      {
         for (int j = M_srcf->ig[i]; j < M_srcf->ig[i + 1]; j++)
         {
            int k = M_srcf->jg[j];
            M_srf->gg[adr[k]] = M_srcf->u[j];
            adr[k]++;
         }
      }
   }

   Matrix* MakeSparseProfileFormatFromRCF(Matrix* M)
   {
      Matrix *N = new Matrix;
      N->di.resize(M->dim);
      N->dim = M->dim;
      for (int i = 0; i < M->dim; i++)
         N->di[i] = M->di[i];

      N->l.clear();
      N->u.clear();
      N->l.reserve(M->ig[M->dim] * 2);
      N->u.reserve(M->ig[M->dim] * 2);
      N->ig.resize(M->dim + 1);

      int col = 0, count = 0;
      for (int i = 0; i < M->dim; i++)
      {
         N->ig[i] = count;
         col = M->jg[M->ig[i]];
         for (int j = M->ig[i]; j < M->ig[i + 1]; j++)
         {
            // i - row, _ja[j] - col
            if (col == M->jg[j])
            {
               N->l.push_back(M->l[j]);
               N->u.push_back(M->u[j]);
               count++;
               col++;
            }
            else
            {
               while (col != M->jg[j])
               {
                  N->l.push_back(0.0);
                  N->u.push_back(0.0);
                  count++;
                  col++;
               }
               N->l.push_back(M->l[j]);
               N->u.push_back(M->u[j]);
               count++;
               col++;
            }
         }
         while (col != i)
         {
            N->l.push_back(0.0);
            N->u.push_back(0.0);
            count++;
            col++;
         }
      }
      N->ig[M->dim] = count;

      N->format = SparseProfile;
      return N;
   }

   void AddElement(Matrix* M, int i, int j, real elem)
   {
      switch (M->format)
      {
      case Dense:
            M->dense[i][j] += elem;
         break;
      case SparseProfile:

         break;
      case SparseRowColumn:
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
         break;

      }

         
   }

   void MatxVec(std::vector<real>& v, Matrix* M, std::vector<real>& b) // v = M*b
   {
      switch (M->format)
      {
      case Dense:
         for (size_t i = 0; i < M->dim; i++)
         {
            real sum = 0;
            for (int j = 0; j < M->dim; j++)
               sum += M->dense[i][j] * b[j];
            v[i] = sum;
         }
         break;
      case SparseProfile:

         break;
      case SparseRowColumn:
         for (int i = 0; i < M->dim; i++)
            v[i] = M->di[i] * b[i];

         for (int i = 0; i < M->dim; i++)
            for (int j = M->ig[i]; j < M->ig[i + 1]; j++) // -1?
            {
               v[i] += M->l[j] * b[M->jg[j]];
               v[M->jg[j]] += M->u[j] * b[i];
            }
         break;
      case SparseRow:
         for (int i = 0; i < M->dim; i++)
         {
            real sum = 0.;
            for (int k = M->ig[i]; k < M->ig[i + 1]; k++)
               sum += M->gg[k] * b[M->jg[k]];

            v[i] = sum;
         }

         break;
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

   void SLAEResidualOutput(std::vector<real>& q, maths::Matrix* M, std::vector<real>& b)
   {

      std::vector<real> y, & x = q;
      real res = 0;
      y.resize(M->dim, 0);
      MatxVec(y, M, x);
      for (size_t i = 0; i < M->dim; i++)
         y[i] -= b[i];
      res = sqrt(scalar(y, y) / scalar(b, b));

      std::cout << "res: " << res << '\n';
      std::cout << "(N = " << M->dim << "): ";
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

         real xsum = 0., rsum = 0.; 
         for (i = 0; i < M->dim; i++)
         {
            xsum += alpha * z[i];
            rsum += alpha * p[i];
         }
         x[i] += xsum;
         r[i] -= rsum;
         MatxVec(Ar, M, r);
         beta = -scalar(p, Ar) / skp;
         for (i = 0; i < M->dim; i++)
         {
            z[i] = r[i] + beta * z[i];
            p[i] = Ar[i] + beta * p[i];
         }
         res = sqrt(scalar(r, r)) / b_norm;
      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void SolveSLAE_PARDISO(Matrix* M, std::vector<real>& q, std::vector<real>& b)
   {
      MKL_INT64 n = M->dim;
      MKL_INT64 mtype = 2; // real and symmetric positive definite
      MKL_INT64 nrhs = 1;
      void* pt[64];
      MKL_INT64 maxfct = 1;
      MKL_INT64 mnum = 1;
      MKL_INT64 msglvl = 1;
      MKL_INT64 phase = 13;
      MKL_INT64* perm = NULL;
      MKL_INT64 iparam[64];
      MKL_INT64 info = -100;

      std::vector<MKL_INT64> ia, ja;

      ia.resize(M->ig.size());
      ja.resize(M->jg.size());
      for (int i = 0; i < M->ig.size(); i++) // MKL_INT64 <- int
         ia[i] = M->ig[i];
      for (int i = 0; i < M->jg.size(); i++)
         ja[i] = M->jg[i];
      for (int i = 0; i < q.size(); i++)
         q[i] = 0.;

      pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n,
          M->gg.data(), ia.data(), ja.data(), perm,
          &nrhs, iparam, &msglvl, b.data(), q.data(), &info);
      //SLAEResidualOutput(q, M, b);

   }

   void SolveSLAE_LU(Matrix*& LU, Matrix* A, std::vector<real>& q, std::vector<real>& b)
   {
      //Matrix* LU = *LUp;

      std::vector<real> y, &x = q;
      x.clear();
      x.resize(A->dim, 0);
      y.resize(A->dim, 0);
      std::clock_t start = clock();

      if (LU == nullptr)
         MakeLUFromRCF(LU, A);

      SolveForL(y, b, LU);
      SolveForU(x, y, LU);
      // L 

      std::clock_t end = clock();

      SLAEResidualOutput(q, A, b);
      std::cout << "Work time (N = " << A->dim << "): " << end - start << "ms [" << (end - start) / 1000 << "s] (" << (end - start) / 60000 << "m)\n";
}

   void MakeLUFromRCF(Matrix*& LU, Matrix* A)
   {
      switch (A->format)
      {
      case Dense:
      {
         LU = MakeDenseFormat(A->dim);
         for (int i = 0; i < LU->dim; i++)
         {
            for (int j = 0; j < LU->dim; j++)
            {
               if (i <= j) // U
               {
                  real sum = 0.;
                  for (int k = 0; k < i; k++)
                     sum += LU->dense[i][k] * LU->dense[k][j];
                  LU->dense[i][j] = A->dense[i][j] - sum;

               }
               else // L
               {
                  real sum = 0.;
                  for (int k = 0; k < j; k++)
                     sum += LU->dense[i][k] * LU->dense[k][j];
                  LU->dense[i][j] = (A->dense[i][j] - sum) / LU->dense[j][j];
               }
            }
         }
         break;
      }
      case SparseProfile:
      case SparseRow: 
         break;

      case SparseRowColumn:
      {
         LU = MakeSparseProfileFormatFromRCF(A);

         for (int i = 1; i < LU->dim; i++)
         {
            int j0 = i - (LU->ig[i + 1] - LU->ig[i]);
            for (int ii = LU->ig[i]; ii < LU->ig[i + 1]; ii++)
            {
               int j = ii - LU->ig[i] + j0;
               double sum_l = 0, sum_u = 0;
               if (LU->ig[j] < LU->ig[j + 1])
               {
                  int j0j = j - (LU->ig[j + 1] - LU->ig[j]);
                  int jjbeg = j0 < j0j ? j0j : j0; // max (j0, j0j)
                  int jjend = j < i - 1 ? j : i - 1; // min (j, i - 1)
                  for (int k = 0; k < jjend - jjbeg; k++)
                  {
                     int ind_prev = LU->ig[j] + jjbeg - j0j + k;
                     int ind_now = LU->ig[i] + jjbeg - j0 + k;
                     sum_l += LU->u[ind_prev] * LU->l[ind_now];
                     sum_u += LU->u[ind_now] * LU->l[ind_prev];
                  }
               }
               LU->l[ii] -= sum_l;
               LU->u[ii] -= sum_u;
               if (abs(LU->di[j]) < 1e-12) // matrix hasn't LU
               {
                  throw std::exception("LU: zero devision!");
               }
               LU->u[ii] /= LU->di[j];
               LU->di[i] -= LU->l[ii] * LU->u[ii];
            }
         }

         break;
      }
      }
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
      SLAEResidualOutput(q, A, b);
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

      for (; abs(res - prevres) < 1e-14 && k < 100000; k++)
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

      SLAEResidualOutput(q, M, b);
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
      switch (M->format)
      {
      case Dense:
         for (int j = 0; j < M->dim; j++)
         {
            b[j] -= b[i] * M->dense[i][j];
            M->dense[i][j] = 0;
         }
         M->dense[i][i] = 1.;
         break;
      case SparseProfile:

         break;
      case SparseRowColumn:
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
         break;

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
      M->format = SparseRowColumn;

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

   void SolveForL(std::vector<real>& q, std::vector<real>& b, Matrix* M) // y = 1/L * b // q = 1/L * b // forward
   {
      switch (M->format)
      {
         case Dense:
            for (int i = 0; i < M->dim; i++)
            {
               real sum = 0;
               for (int j = 0; j < i; j++)
                  sum += q[j] * M->dense[i][j];
               q[i] = b[i] - sum;
            }
            break;
         case SparseProfile:
            for (int i = 0; i < M->dim; i++)
            {
               int k = i - (M->ig[i + 1] - M->ig[i]);
               real sum = 0;
               for (int j = M->ig[i]; j < M->ig[i + 1]; j++, k++)
                  sum += q[k] * M->l[j];
                //
               q[i] = (b[i] - sum) / M->di[i];
            }
            break;
         case SparseRowColumn:
            for (int k = 1, k1 = 0; k < M->dim; k++, k1++)
            {
               double sum = 0;
               for (int i = M->ig[k1]; i < M->ig[k]; i++)
                  sum += M->l[i] * q[M->jg[i]];

               q[k1] = b[k1] - sum;
            }
            break;

      }
   }
   void SolveForU(std::vector<real>& q, std::vector<real>& b, Matrix* M) // x = 1/U * y // q = 1/U * b // back
   {
      switch (M->format)
      {
      case Dense:
         for (int i = M->dim - 1; i >= 0; i--)
         {
            real sum = b[i];
            for (int j = i + 1; j < M->dim; j++)
               sum -= q[j] * M->dense[i][j];
            q[i] = sum / M->dense[i][i];
         }
         break;

      case SparseProfile:
         for (int i = M->dim - 1; i >= 0; i--)
         {
            //real sum = b[i];
            int j = i - (M->ig[i + 1] - M->ig[i]);
            q[i] += b[i];
            for (int k = M->ig[i]; k < M->ig[i + 1]; k++, j++)
               q[j] -= M->u[k] * q[i];
         }

         break;

      case SparseRowColumn: // неправильно
         for (int i = M->dim - 1; i >= 0; i--)
         {
            real sum = b[i];

            for (int j = M->ig[i]; j < M->ig[i + 1]; j++) // -1?
               sum -= M->u[j] * q[M->jg[j]];

            q[i] = sum / M->di[i];
         }
         break;

      }
   }
}