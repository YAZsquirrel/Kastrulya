#include "FEM.h"
#include <string>
FEM::FEM()
{
   std::ifstream ftime("TimeGridDescr.txt");
   ftime >> t_last >> dt >> tr >> u0 >> utest;
   ftime.close();
   dt = dt > 1e-10 ? t_last / dt : 0.0;
   mesh = new Mesh();
   mesh->MakeMesh();
   num_of_knots = mesh->knots.size();
   num_of_FE = mesh->elems.size();
#ifndef DENSE
   A = MakeSparseRowColumnFormat(4, mesh->knots.size(), mesh);
    //A => G, M, Gx
   {
      G = new Matrix();
      Gx = new Matrix();
      M = new Matrix();
      Gx->dim = G->dim = M->dim = A->dim;
      M->ig.resize(A->ig.size());
      G->ig.resize(A->ig.size());
      Gx->ig.resize(A->ig.size());
      M->jg.resize(A->jg.size());
      G->jg.resize(A->jg.size());
      Gx->jg.resize(A->jg.size());
      copy(M->ig, A->ig);
      copy(G->ig, A->ig);
      copy(Gx->ig, A->ig);
      copy(M->jg, A->jg);
      copy(G->jg, A->jg);
      copy(Gx->jg, A->jg);
      G->di.resize(A->dim);
      Gx->di.resize(A->dim);
      M->di.resize(A->dim);
      M->l.resize(A->l.size());
      G->l.resize(A->l.size());
      Gx->l.resize(A->l.size());
      M->u.resize(A->u.size());
      G->u.resize(A->u.size());
      Gx->u.resize(A->u.size());
      M->format = G->format = Gx->format = SparseRowColumn;
   }
#else
   A = MakeDenseFormat(mesh->knots.size());
   G = MakeDenseFormat(mesh->knots.size());
   Gx = MakeDenseFormat(mesh->knots.size());
   M = MakeDenseFormat(mesh->knots.size());

#endif // DEBUG


   q1.resize(num_of_knots, 0.);
   q2.resize(num_of_knots, 0.);
   b.resize(num_of_knots, 0.);
   d.resize(num_of_knots, 0.);
}

void FEM::SolveElliptic()
{
   CreateSLAE(false);
   SolveSLAE_LOS(A, q1, b);
}

void FEM::SolveParabolic()
{
   //CreateSLAE();
   if (utest > 0)
#ifdef DEBUG2
   {
      for (int i = 0; i < q1.size(); i++)
         q1[i] = bound1func(mesh->knots[i], 0, utest);
      CreateSLAE(false);
   }
#else
      SolveElliptic();
#endif // DEBUG
   else {
      for (int i = 0; i < q1.size(); i++)
         q1[i] = u0;
      CreateSLAE(false);
   }
   int tn = 0;
   std::string str = "./Results/Result_layer_";

   std::vector<real> zero;
   zero.resize(num_of_knots, 0.0);

   real t = 0;
   for (; abs(t - dt - t_last) > 1e-12; t += dt, tn++)
   {
      copy(b, zero);
      if (tn > 0)
      {
         for (int i = 0; i < num_of_FE; i++)
         {
            CreateM(mesh->elems[i]);
            Createb(mesh->elems[i], t);
         }
         AssembleMatricies(true, t);
         //SolveSLAE_LOSnKholessky(A, q1, d);
         //SolveSLAE_Relax(A, q1, d, 1.6);
#ifdef DENSE
         SolveSLAE_LU(A, q1, d);
         //SolveSLAE_LOS(A, q1, d);
#else
         SolveSLAE_LU(A, q1, d);
         //SolveSLAE_LOS(A, q1, d);
#endif // DENSE

         //WriteMatrix(A);
         //copy(q1, q2);

      }  
      std::cout << "Layer " << tn << ", Time: " << t << "s, done." << std::endl;
      std::ofstream out(str + std::to_string(tn) + ".txt", std::ofstream::in);
      out.close();
      out.open(str + std::to_string(tn) + ".txt", std::ofstream::trunc);
      Output(out);
      out.close();
   }

   if (tn == 0)
   {
      std::ofstream out(str + std::to_string(tn) + ".txt", std::ofstream::in);
      out.close();
      out.open(str + std::to_string(tn) + ".txt", std::ofstream::trunc);
      Output(out);
      out.close();
   }
}

void FEM::Output(std::ofstream& out)
{
   out.setf(std::ios::right);
   out.width(15);
   out << "\nx" << std::fixed;
   out.width(15);
   out << "y";
   out.width(15);
   //out << "z";
   //out.width(15);
   out << "q";
   out.width(15);
   out << "\n";
   out << std::setprecision(7);

   for (int i = 0; i < num_of_knots; i++)
   {
      out << std::defaultfloat;
      out << mesh->knots[i].x;
      out.width(15);
      out << mesh->knots[i].y;
      out.width(30);
      //out << mesh->knots[i].z;
      //out.width(15);
      out << std::scientific;
      out.width(15);
      out << q1[i];
      out.width(15);
      out << "\n";
   }

   std::ofstream out2("speed.txt", std::ofstream::in);
   out.close();
   out.open("speed.txt", std::ofstream::trunc);

   out.setf(std::ios::right);
   out.width(15);
   out << "\nx" << std::fixed;
   out.width(15);
   out << "y";
   out.width(15);
   //out << "z";
   //out.width(15);
   out << "v";
   out.width(15);
   out << "\n";
   out << std::setprecision(7);

   for (int i = 0; i < num_of_knots; i++)
   {
      knot v = mesh->GetVelocity(mesh->knots[i], 1);
      out << std::defaultfloat;
      out << mesh->knots[i].x;
      out.width(15);
      out << mesh->knots[i].y;
      out.width(30);
      //out << mesh->knots[i].z;
      //out.width(15);
      out << std::scientific;
      out.width(15);
      out << sqrt(v.x * v.x + v.y * v.y);
      out.width(15);
      out << "\n";
   }


}

void FEM::AddFirstBounds(real time)
{
   for (auto& cond : mesh->bounds1)
   {
      for (int i = 0; i < 2; i++)
      {
         switch (A->format)
         {
            case Dense:
               for (size_t ii = 0; ii < A->dim; ii++)
                  A->dense[i][ii] = 0.;
               A->dense[i][i] = 1.;

               break;
            case SparseProfile:
               break;

            case SparseRowColumn:
               A->di[cond.knots_num[i]] = 1.;
               for (int j = A->ig[cond.knots_num[i]]; j < A->ig[cond.knots_num[i] + 1]; j++)
                  A->l[j] = 0.;
               for (int ii = 0; ii < A->dim; ii++)                // ???? ?? ????????
                  for (int j = A->ig[ii]; j < A->ig[ii + 1]; j++)   // ???? ????????? ? ???????
                     if (A->jg[j] == cond.knots_num[i])          // ? ?????? ?????? ????????
                        A->u[j] = 0.;
               break;
         }

#ifdef DEBUG2
         if (time > 1e-10)
            d[cond.knots_num[i]] = bound1func(mesh->knots[cond.knots_num[i]], time, cond.n_test);
         else
            b[cond.knots_num[i]] = bound1func(mesh->knots[cond.knots_num[i]], time, cond.n_test);
#else
         if (time > 1e-10)
            d[cond.knots_num[i]] = cond.value1;
         else
            b[cond.knots_num[i]] = cond.value1;
#endif // DEBUG

         if (time > 1e-10)
            MatSymmetrisation(A, d, cond.knots_num[i]);
         else 
            MatSymmetrisation(A, b, cond.knots_num[i]);
      }
   }
}

void FEM::AddSecondBounds(real time)
{
   for (auto& bound : mesh->bounds2)
   {
      real r1 = std::min(mesh->knots[bound.knots_num[0]].x, mesh->knots[bound.knots_num[1]].x),
           r2 = std::max(mesh->knots[bound.knots_num[0]].x, mesh->knots[bound.knots_num[1]].x),
           z1 = std::min(mesh->knots[bound.knots_num[0]].y, mesh->knots[bound.knots_num[1]].y),
           z2 = std::max(mesh->knots[bound.knots_num[0]].y, mesh->knots[bound.knots_num[1]].y),
           h, h2;

      bool isRadiusAxis = r2 - r1 > 1e-14;

      if (isRadiusAxis)
      {
         h = r2 - r1;
         h2 = pow(h, 2.);
         real h2_4 = h*h/4.,
              hr_3 = h*r1/3.;
         localMb[0][0] = hr_3      + h2_4 / 3.;  localMb[0][1] = hr_3 / 2. + h2_4 / 3.;
         localMb[1][0] = hr_3 / 2. + h2_4 / 3.;  localMb[1][1] = hr_3      + h2_4;
      }
      else
      {
         h = z2 - z1;
         h2 = pow(h, 2.);
         real rh_3 = r1 * h / 3.;
         localMb[0][0] = rh_3;       localMb[0][1] = rh_3 / 2.;
         localMb[1][0] = rh_3 / 2.;  localMb[1][1] = rh_3;
      }
      //real ratio = len;
      
#ifdef DEBUG2
      for (int i = 0; i < 2; i++)
      {
         if (time > 1e-10)
            d[bound.knots_num[i]] += bound2func(mesh->knots[bound.knots_num[i]], time, bound.n_test) * (localMb[i][0] + localMb[i][1]);
         else
            b[bound.knots_num[i]] += bound2func(mesh->knots[bound.knots_num[i]], time, bound.n_test) * (localMb[i][0] + localMb[i][1]);
      }
#else
      for (int i = 0; i < 2; i++)
      {
         if (time > 1e-10)
            d[bound.knots_num[i]] += bound.value1 * (localMb[i][0] + localMb[i][1]);
         else
            b[bound.knots_num[i]] += bound.value1 * (localMb[i][0] + localMb[i][1]);
      }
#endif // DEBUG
   }
}

void FEM::AddThirdBounds(real time)
{
   for (auto& bound : mesh->bounds3)
   {
    real r1 = std::min(mesh->knots[bound.knots_num[0]].x, mesh->knots[bound.knots_num[1]].x),
         r2 = std::max(mesh->knots[bound.knots_num[0]].x, mesh->knots[bound.knots_num[1]].x),
         z1 = std::min(mesh->knots[bound.knots_num[0]].y, mesh->knots[bound.knots_num[1]].y),
         z2 = std::max(mesh->knots[bound.knots_num[0]].y, mesh->knots[bound.knots_num[1]].y),
         h, h2;

      bool isRadiusAxis = r2 - r1 > 1e-10;

      if (isRadiusAxis)
      {
         h = r2 - r1;
         h2 = pow(h, 2.);
         real h2_4 = h*h/4.,
              hr_3 = h*r1/3.;
         localMb[0][0] = hr_3      + h2_4 / 3.;  localMb[0][1] = hr_3 / 2. + h2_4 / 3.;
         localMb[1][0] = hr_3 / 2. + h2_4 / 3.;  localMb[1][1] = hr_3      + h2_4;
      }
      else
      {
         h = z2 - z1;
         h2 = pow(h, 2.);
         real rh_3 = r1 * h / 3.;
         localMb[0][0] = rh_3;       localMb[0][1] = rh_3 / 2.;
         localMb[1][0] = rh_3 / 2.;  localMb[1][1] = rh_3;
      }
      //real ratio = len;

#ifdef DEBUG2
      for (int i = 0; i < 2; i++)
      {
         if (time > 1e-10)
            d[bound.knots_num[i]] += bound3func(mesh->knots[bound.knots_num[i]], time, bound.n_test)
                                    * bound3funcbeta(mesh->knots[bound.knots_num[i]], time, bound.n_test)
                                    * (localMb[i][0] + localMb[i][1]);
         else
            b[bound.knots_num[i]] += bound3func(mesh->knots[bound.knots_num[i]], time, bound.n_test)
                                    * bound3funcbeta(mesh->knots[bound.knots_num[i]], time, bound.n_test)
                                    * (localMb[i][0] + localMb[i][1]);

         for (int j = 0; j < 2; j++)
            AddElement(A, bound.knots_num[i], bound.knots_num[j], bound3funcbeta(mesh->knots[bound.knots_num[i]], time, bound.n_test) * localMb[i][j]);
      }
#else
      for (int i = 0; i < 2; i++)
      {
         if (time > 1e-10)
            d[bound.knots_num[i]] += bound.value1 * bound.value2 * (localMb[i][0] + localMb[i][1]);
         else
            b[bound.knots_num[i]] += bound.value1 * bound.value2 * (localMb[i][0] + localMb[i][1]);

         for (int j = 0; j < 2; j++)
            AddElement(A, bound.knots_num[i], bound.knots_num[j], bound.value2 * localMb[i][j]);
      }
#endif // DEBUG

   }
}

void FEM::CreateSLAE(bool isTimed)
{
   element elem;
   for (int i = 0; i < num_of_FE; i++)
   {
      elem = mesh->elems[i];
      CreateG(elem);
      CreateM(elem);
      if (elem.n_mat == 0)
         CreateExtraG(elem);
      Createb(elem, 0);
      AddToGlobalMatricies(elem);
   }

   AssembleMatricies(isTimed, 0);
}

void FEM::AssembleMatricies(bool isTimed, real time)
{
   real mulM = isTimed ? 1. / dt : 1.;
   real mulGx = isTimed ? 1.0 : 0.0;
   switch (A->format)
   {
   case Dense:
      for (size_t i = 0; i < A->dim; i++)
         for (size_t j = 0; j < A->dim; j++)
            A->dense[i][j] = mulM * M->dense[i][j] + G->dense[i][j] + mulGx * Gx->dense[i][j];

      break;
   case SparseProfile:
      break;

   case SparseRowColumn:
      for (int i = 0; i < A->l.size(); i++)
      {
         A->l[i] = mulM * M->l[i] + G->l[i] + mulGx * Gx->l[i];
         A->u[i] = mulM * M->u[i] + G->u[i] + mulGx * Gx->u[i];
      }
      for (int i = 0; i < A->dim; i++)
         A->di[i] = mulM * M->di[i] + G->di[i] + mulGx * Gx->di[i];
      break;
   }


   if (isTimed)
   {
      std::vector<real> Mq;
      Mq.resize(num_of_knots);
      MatxVec(Mq, M, q1);
      for (int i = 0; i < Mq.size(); i++)
         d[i] = b[i] + mulM * Mq[i];
   }
   AddSecondBounds(time);
   AddThirdBounds(time);
   AddFirstBounds(time);
}

void FEM::AddToGlobalMatricies(element& elem)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
      {
         //localA[i][j] = localG[i][j] + localM[i][j];
         AddElement(M, elem.knots_num[i], elem.knots_num[j], localM[i][j]);
         AddElement(G, elem.knots_num[i], elem.knots_num[j], localG[i][j]);
         AddElement(Gx, elem.knots_num[i], elem.knots_num[j], localGx[i][j]);
         //AddElement(A, elem.knots_num[i], elem.knots_num[j], localA[i][j]);
      }

//AddLocal(A, elem->knots_num, localA, 1);
}

void FEM::CreateM(element& elem)
{
   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real h = r2 - r1, t = z2 - z1;
   real h2t_72 = h * h * t / 72., 
        htr_36 = h * t * r1 / 36.;

   real localMz[4][4] = {{4,2,2,1},
                         {2,4,1,2},
                         {2,1,4,2},
                         {1,2,2,4}},
        localMr[4][4] = {{2,2,1,1},
                         {2,6,1,3},
                         {1,1,2,2},
                         {1,3,2,6} };
         

   //localC[0][0] = 
   //localC[1][1] = h2t_72 * 2. + htr_36 * 4.;       // +
   //localC[2][2] = 
   //localC[3][3] = h2t_72 * 6. + htr_36 * 4.;       // +
   //localC[0][3] = localC[3][0] =                            // +
   //localC[1][2] = localC[2][1] = h2t_72      + htr_36;  // +
   //localC[0][2] = localC[2][0] =                            //
   //localC[1][3] = localC[3][1] = h2t_72 * 2. + htr_36 * 2.;  //
   //localC[0][1] = localC[1][0] = h2t_72      + htr_36 * 2.;  // +
   //localC[2][3] = localC[3][2] = h2t_72 * 3. + htr_36 * 2.;  // +

   // 2 1 2 1
   // 1 2 1 2
   // 2 1 6 3
   // 1 2 3 6

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
      {
         localC[i][j] = h2t_72 * localMr[i][j] + htr_36 * localMz[i][j];
         localM[i][j] = elem.gam * localC[i][j];
      }        
}

void FEM::CreateG(element& elem)
{
   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real hr = r2 - r1, hz = z2 - z1;

   real h = hz * (2 * r1 + hr) / hr / 12.,
        p = hr * r1 / hz / 6.,
        q = hr * hr / hz / 12.;
   real localH[4][4]{ {2,-2,1,-1},
                      {-2,2,-1,1},
                      {1,-1,2,-2},
                      {-1,1,-2,2} },

        localP[4][4]{ {2,1,-2,-1},
                      {1,2,-1,-2},
                      {-2,-1,2,1},
                      {-1,-2,1,2} },

        localQ[4][4]{ {1,1,-1,-1},
                      {1,3,-1,-3},
                      {-1,-1,1,1},
                      {-1,-3,1,3} };

   //localG[0][0] = 
   //localG[2][2] = 2. * h + 2. * p +      q;
   //localG[1][1] = 
   //localG[3][3] = 2. * h + 2. * p + 3. * q;
   //localG[1][2] = localG[2][1] = 
   //localG[0][3] = localG[3][0] = - h - p - q;
   //localG[0][1] = localG[1][0] =
   //localG[2][3] = localG[3][2] =  -2 * h + p + q;
   //localG[1][3] = localG[3][1] =  h - 2. * p - 3. * q;
   //localG[0][2] = localG[2][0] =  h - 2. * p - q;



   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localG[i][j] = elem.lam * (h * localH[i][j] + p * localP[i][j] + q * localQ[i][j]);//* Integrate(Gij, i, j, elem->knots_num);

}

void FEM::CreateExtraG(element& elem)
{
   knot k1 =  mesh->knots[elem.knots_num[0]];
   knot k2 =  mesh->knots[elem.knots_num[1]];
   knot k3 =  mesh->knots[elem.knots_num[2]];
   knot k4 =  mesh->knots[elem.knots_num[3]];

   knot v4[4] = {mesh->GetVelocity(k1, 1),
                mesh->GetVelocity(k2, 1),
                mesh->GetVelocity(k3, 1),
                mesh->GetVelocity(k4, 1)};

   knot c = mesh->intersection(k1, k4, k2, k3);
   knot v = mesh->GetVelocity(c, 1);

   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real hr = r2 - r1, hz = z2 - z1;

   real h = hr * hz / 36.,
        p = hz * r1 / 12.,
        q = hr * r1 / 12.,
        t = hr * hr / 24.;
   // my
   real localH[4][4]{{-2,2,-1,1}, 
                     {-4,4,-2,2},
                     {-1,1,-2,2},
                     {-2,2,-4,4}},
        localP[4][4]{{-2,2,-1,1},
                     {-2,2,-1,1},
                     {-1,1,-2,2},
                     {-1,1,-2,2} },
        localQ[4][4]{{-2,-1,2,1},
                     {-1,-2,1,2},
                     {-2,-1,2,1},
                     {-1,-2,1,2} },
        localT[4][4]{{-1,-1,1,1},
                     {-1,-3,1,3},
                     {-1,-1,1,1},
                     {-1,-3,1,3} };


   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         //localGx[i][j] = elem.gam * (v4[i].x * h * localH[i][j] + v4[i].x * p * localP[i][j] + v4[i].y * q * localQ[i][j] + v4[i].y * t * localT[i][j]);
         localGx[i][j] = elem.gam * (v.x * h * localH[i][j] + v.x * p * localP[i][j] + v.y * q * localQ[i][j] + v.y * t * localT[i][j]);
}

void FEM::Createb(element& elem, real time)
{
   real localb[4]{};
   real f_[4]{};


#ifdef DEBUG2
   for (int i = 0; i < 4; i++)
      f_[i] = f(mesh->knots[elem.knots_num[i]], elem, time);
#else
   for (int i = 0; i < 4; i++)
      f_[i] = 0;//elem.f;
#endif // DEBUG

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localb[i] += localC[i][j] * f_[j];

   for (int i = 0; i < 4; i++)
      b[elem.knots_num[i]] += localb[i];
}
