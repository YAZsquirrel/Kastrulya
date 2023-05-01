#include "Math_structs.h"
#include <set>
#include <array>
#include <cmath>
#include <atomic>
#include <numeric>
#include <array>
#include <mkl_pardiso.h>
#include <unordered_set>

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

      Matrix* M = new Matrix(MatrixFormat::SparseRowColumn);
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

   void MakeSparseRowFormatFromRCF(Matrix* M_rcf, Matrix*& M_srf)
   {
       M_srf = new Matrix(MatrixFormat::SparseRow);

       int n = M_srf->dim = M_rcf->dim;

       M_srf->ig.resize(n + 1);
       M_srf->jg.resize(M_rcf->ig[n] + n);
       M_srf->gg.resize(M_rcf->ig[n] + n, 0);
       copy(M_srf->di, M_rcf->di);

       ConvertFromSRCFToSRF(M_rcf, M_srf);
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

      for (int i = 0; i < N; i++)
         adr[i] = M_srf->ig[i];

      int jbuf, gbuf;
      for (int i = 0; i < N; i++)
      {
         M_srf->jg[adr[i]] = i;
         adr[i]++;
         for (int jrc = M_srcf->ig[i], jr = M_srf->ig[i]; jrc < M_srcf->ig[i + 1]; jrc++, jr++)
         {
            jbuf = M_srcf->jg[jrc];
            M_srf->jg[adr[jbuf]] = i;
            adr[jbuf]++;
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
         adr[i] = M_srf->ig[i];

      int jbuf;
      double gbuf;
      for (int i = 0; i < N; i++)
      {
         M_srf->gg[adr[i]] = M_srcf->di[i];
         adr[i]++;
         for (int jrc = M_srcf->ig[i], jr = M_srf->ig[i]; jrc < M_srcf->ig[i + 1]; jrc++, jr++)
         {

            jbuf = M_srcf->jg[jrc];
            gbuf = M_srcf->l[jrc];
            M_srf->gg[adr[jbuf]] = gbuf;
            adr[jbuf]++;
         }
      }
   }

   Matrix* MakeSparseProfileFormatFromRCF(Matrix* M)
   {
      Matrix *N = new Matrix(MatrixFormat::SparseProfile);
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
            // i - row, _ja[ii] - col
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

   void MatxVec(std::vector<real>& v, Matrix* M, std::vector<real>& b) // v = M_rcf*b
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
            for (int ii = M->ig[i]; ii < M->ig[i + 1]; ii++) // -1?
            {
               v[i] += M->l[ii] * b[M->jg[ii]];
               v[M->jg[ii]] += M->u[ii] * b[i];
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

   real SLAEResidualOutput(std::vector<real>& q, maths::Matrix* M, std::vector<real>& b)
   {
      std::vector<real> y, & x = q;
      y.resize(M->dim, 0);
      MatxVec(y, M, x);
      for (size_t i = 0; i < M->dim; i++)
         y[i] -= b[i];
      return sqrt(scalar(y, y)) / sqrt(scalar(b, b));

      //std::cout << "res: " << res << '\n';
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
      }
      std::cout << "iter: " << k << " Residual: " << res << std::endl;
   }

   void SolveSLAE_PARDISO(Matrix* M, std::vector<real>& q, std::vector<real>& b)
   {

      std::ofstream logfile;
      logfile.open("pardiso64.log");
      if (!logfile)
      {
         std::cerr << "Cannot open pardiso64.log\n";
         logfile.close();
      }

      MKL_INT64 n = M->dim;
      MKL_INT64 mtype = 2; // real and symmetric positive definite
      MKL_INT64 nrhs = 1;
      void* pt[64]{};
      MKL_INT64 maxfct = 1;
      MKL_INT64 mnum = 1;
      MKL_INT64 msglvl = 0;
      MKL_INT64 phase = 13;
      MKL_INT64* perm = new MKL_INT64[n]{};
      MKL_INT64 iparam[64]{};
      MKL_INT64 info = -100;

      std::vector<MKL_INT64> ia, ja;

      ia.resize(M->ig.size());
      ja.resize(M->jg.size());
      for (int i = 0; i < M->ig.size(); i++) // MKL_INT64 <- int
         ia[i] = M->ig[i] + 1;
      for (int i = 0; i < M->jg.size(); i++)
         ja[i] = M->jg[i] + 1;
      //for (int i = 0; i < q.size(); i++)
      //   q[i] = 0.;

      pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n,
         M->gg.data(), ia.data(), ja.data(),
                  perm, &nrhs, iparam, &msglvl, 
                  (void*)b.data(), (void*)q.data(),
                  &info);
                  
      logfile << info << '\n';
      logfile.close();
      logfile.clear();
      if (perm) { delete[] perm; perm = NULL; }
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
                  int jjend = j < i - 1 ? j : i - 1; // min (ii, i - 1)
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

   void SolveSLAE_predet_LOS(Matrix* A, std::vector<real>& q, std::vector<real>& b)
   {
      std::vector<real> z, r, p, AQr, t, SAQr, &x = q;
      z.resize(A->dim);
      r.resize(A->dim);
      p.resize(A->dim);
      AQr.resize(A->dim);
      SAQr.resize(A->dim);
      t.resize(A->dim);
      real res, alpha, beta, pp, eps = 1e-14;
      int i, k;
      //x = q;

      Matrix* SQ = MakeHolessky(A);

      MatxVec(t, A, x);  // A*x0                          
      real b_norm;

      for (int i = 0; i < A->dim; i++)    //b - A*x0     
         t[i] = b[i] - t[i];

      SolveForL(r, t, SQ); // r0 = 1/S * (b - Ax0) <-> S*r0 = b - A*x0 
      SolveForU(z, r, SQ);  // z0 = 1/Q * r0 <-> Q*z0 = r0

      MatxVec(t, A, z);    // A*z0                         
      SolveForL(p, t, SQ); // p0 = 1/S * A*z0 <-> S*p0 = A*z0

      b_norm = sqrt(scalar(b, b));
      real res0 = res = scalar(r, r); //sqrt(scalar(r, r)) / b_norm;

      for (k = 1; k < 100000 && res / res0 > eps * eps; k++)
      {
         for (int i = 0; i < A->dim; i++)
            SAQr[i] = AQr[i] = t[i] = 0.;

         pp = scalar(p, p);
         alpha = scalar(p, r) / pp;          // a = (pk0, rk0) / (pk0, pk0)

         for (int i = 0; i < A->dim; i++)     // xk1 = xk0 + a*zk0
         {                                    // rk1 = rk0 - a*pk0 
            x[i] += alpha * z[i];               
            r[i] -= alpha * p[i];              
         }
                        
         SolveForU(t, r, SQ);                 // t = 1/Q * rk <-> Q*t = rk
         MatxVec(AQr, A, t);                  // A * t1 = AQr -> AQr = A * 1/Q * rk
         SolveForL(SAQr, AQr, SQ);             // SAQr = 1/S AQr <-> S*SAQr = AQr -> SAQr = 1/S * A * 1/Q * rk

         beta = -scalar(SAQr, p) / pp;          // beta = (pk0, t) / (pk0, pk0)          

         //SolveForU(t2, r, SQ);                 // Q*zk1 = rk1
         for (int i = 0; i < A->dim; i++)
         {
            z[i] = t[i] + beta * z[i];        // zk1 += beta
            p[i] = SAQr[i] + beta * p[i];
         }
         
         res = scalar(r, r);//sqrt(scalar(r, r)) / b_norm;
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

   Matrix* MakeHolessky(Matrix *A) 
   {

      Matrix* M = new Matrix();
      M->dim = A->dim;
      M->ig.resize(A->ig.size());
      M->jg.resize(A->jg.size());
      M->l.resize(A->l.size());
      M->u.resize(A->u.size());
      M->di.resize(A->di.size());
      M->format = SparseRowColumn;

      copy(M->ig, A->ig);
      copy(M->jg, A->jg);

      for (int i = 0; i < M->dim; i++)
      {
         // Sij, Qji
         for (int ii = M->ig[i]; ii < M->ig[i + 1]; ii++) // ii -> index on row i
         {
            int j = M->jg[ii];

            real sumS = 0.;
            real sumQ = 0.;
            std::vector<int*> ks;
            int is = M->ig[i], ie = M->ig[i + 1];
            int js = M->ig[j], je = M->ig[j + 1];

            for (int ki = is; ki < ie; ki++)
               for (int kj = js; kj < je; kj++)
                  if (M->jg[ki] == M->jg[kj])
                     ks.push_back(new int[2] { ki, kj }); // find intersection of nonzero elem's ii

            for (auto k : ks)
            {
               int ki = k[0], kj = k[1];
               sumS += M->l[ki] * M->u[kj];
               sumQ += M->l[kj] * M->u[ki];
               delete[] k;
            }

            M->l[ii] = (A->l[ii] - sumS) / M->di[j];
            M->u[ii] = (A->u[ii] - sumQ) / M->di[j];

         }
         //Sii
         {
            real sumd = 0.;
            for (int ii = M->ig[i]; ii < M->ig[i + 1]; ii++)
               sumd += M->l[ii] * M->u[ii];

            M->di[i] = sqrt(A->di[i] - sumd);
         }
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
            for (size_t k = 0; k < M->dim; k++)
            {
               real sum = 0.;
               for (size_t i = M->ig[k]; i < M->ig[k + 1]; i++)
                  sum += M->l[i] * q[M->jg[i]];

               q[k] = (b[k] - sum) / M->di[k];
            }

            break;

      }
   }
   void SolveForU(std::vector<real>& q, std::vector<real>& b, Matrix* M) // x = 1/U * y // q = 1/U * b // back
   {
      switch (M->format)
      {
      case Dense:
         for (int i = M->dim - 1; i > 0; i--)
         {
            real sum = b[i];
            for (int j = i + 1; j < M->dim; j++)
               sum -= q[j] * M->dense[i][j];
            q[i] = sum / M->dense[i][i];
         }
         break;

      case SparseProfile:
         for (int i = M->dim - 1; i > 0; i--)
         {
            //real sum = b[i];
            int j = i - (M->ig[i + 1] - M->ig[i]);
            q[i] += b[i];
            for (int k = M->ig[i]; k < M->ig[i + 1]; k++, j++)
               q[j] -= M->u[k] * q[i];
         }

         break;

      case SparseRowColumn: 
         copy(q, b);
         for (int k = M->dim - 1; k >= 0; k--)
         {
            int ii = k;
            real qii = q[k] / M->di[k];
            q[k] = qii;

            for (int i = M->ig[k]; i < M->ig[k + 1]; i++)
               q[M->jg[i]] -= M->u[i] * qii;

         }
         break;

      }
   }
}