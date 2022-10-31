#include "FEM.h"
#include <string>
FEM::FEM()
{
   //std::ifstream fknots("Knots.txt");
   //std::ifstream felems("elements.txt");
   //std::ifstream fbounds1("FirstBounds.txt");
   //std::ifstream fbounds2("SecondBounds.txt");
   //std::ifstream fbounds3("ThirdBounds.txt");
   std::ifstream fparams("EnvParams.txt");
   fparams >> un;
   fparams.close();
   std::ifstream ftime("TimeGridDescr.txt");
   ftime >> t_last >> th >> tr;
   ftime.close();
   th = t_last / th;
   mesh = new Mesh();
   mesh->MakeMesh();
   num_of_knots = mesh->knots.size();
   num_of_FE = mesh->elems.size();
   u0 = mesh->bounds3[0].value1;

   A = MakeSparseFormat(4, mesh->knots.size(), mesh);
   // A => G, M, Gx
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
   }

   q1.resize(num_of_knots, 0.);
   q2.resize(num_of_knots, 0.);
   b.resize(num_of_knots, 0.);
   d.resize(num_of_knots, 0.);

   //Mij = [this](real ksi, real etta, real theta, int i, int j, int knot_num[8])
   //{
   //   for (int ip = 0; ip < 3; ip++)
   //      for (int jp = 0; jp < 3; jp++)
   //         J[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta, theta);
   //
   //   return phi(i, ksi, etta, theta) * phi(j, ksi, etta, theta) * abs(det_J());
   //};
   //
   //Gij = [this](real ksi, real etta, real theta, int i, int j, int knot_num[8])
   //{
   //   for (int ip = 0; ip < 3; ip++)
   //      Jgrad_i[ip] = Jgrad_j[ip] = 0.;
   //
   //   for (int ip = 0; ip < 3; ip++)                                       // | dx/d(ksi)   | dy/d(ksi)   | dz/d(ksi)   |
   //      for (int jp = 0; jp < 3; jp++)                                    // | dx/d(etta)  | dy/d(etta)  | dz/d(etta)  |
   //         J[ip][jp] = prime_by_var(ip, jp, knot_num, ksi, etta, theta);  // | dx/d(theta) | dy/d(theta) | dz/d(theta) |
   //
   //   // J^-1
   //   for (int ip = 0; ip < 3; ip++)
   //      for (int jp = 0; jp < 3; jp++)
   //      {
   //         real min[4]{};
   //         int k = 0;
   //         for (int im = 0; im < 3; im++)
   //            for (int jm = 0; jm < 3; jm++)
   //            {
   //               if (im != ip && jm != jp)
   //                  min[k++] = J[im][jm];
   //            }
   //         reversed_J[jp][ip] = ((ip + jp + 2) % 2 ? -1 : 1) * (min[0] * min[3] - min[1] * min[2]);
   //      }
   //
   //   // grad(phi(ksi, etta, theta))
   //   calc_grad(1, i, ksi, etta, theta);
   //   calc_grad(2, j, ksi, etta, theta);
   //
   //   // J^-1 * grad(phi)
   //   for (int ip = 0; ip < 3; ip++)
   //      for (int jp = 0; jp < 3; jp++)
   //      {
   //         Jgrad_i[ip] += reversed_J[ip][jp] * gradi[jp];
   //         Jgrad_j[ip] += reversed_J[ip][jp] * gradj[jp];
   //      }
   //
   //   // Jgrad_i^T * Jgrad_j
   //   real res = 0;
   //   for (int ip = 0; ip < 3; ip++)
   //      res += Jgrad_i[ip] * Jgrad_j[ip];
   //   return res / abs(det_J());
   //};

}

void FEM::SolveElliptic()
{
   CreateSLAE(false);
   SolveSLAE(A, q1, b);

#ifdef DEBUG
   check_test();
#endif // DEBUG
}

void FEM::SolveParabolic()
{
   //CreateSLAE();
   SolveElliptic();
   int tn = 0;
   std::string str = "./Results/Result_layer_";
   for (real t = 0; t < t_last; t += th, tn++)
   {
      if (tn > 1)
      {
         AssembleMatricies(true);
         SolveSLAE(A, q2, d);
      }
       
      std::ofstream out(str + std::to_string(tn) + ".txt", std::ofstream::in);
      out.close();
      out.open(str + std::to_string(tn) + ".txt", std::ofstream::trunc);
      Output(out);
      out.close();

      copy(q1, q2);
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
}

void FEM::AddFirstBounds()
{
   for (auto& cond : mesh->bounds1)
   {
      for (int i = 0; i < 4; i++)
      {
         A->di[cond.knots_num[i]] = 1.;
         for (int j = A->ig[cond.knots_num[i]]; j < A->ig[cond.knots_num[i] + 1]; j++)
            A->l[j] = 0.;
         for (int ii = 0; ii < A->dim; ii++)                // идем по столбцам
            for (int j = A->ig[ii]; j < A->ig[ii + 1]; j++)   // идем элементам в столбце
               if (A->jg[j] == cond.knots_num[i])          // в нужной строке элемент?
                  A->u[j] = 0.;
         b[cond.knots_num[i]] = cond.value1;

         MatSymmetrisation(A, b, cond.knots_num[i]);
      }
   }
#ifdef DEBUG1
   for (auto cond : mesh->bounds2)
   {
      for (int i = 0; i < 4; i++)
      {
         //A->di[cond.knots_num[i]] = 1e10;
         A->di[cond.knots_num[i]] = 1.;//e10;
         for (int j = A->ig[cond.knots_num[i]]; j < A->ig[cond.knots_num[i] + 1]; j++)
            A->l[j] = 0.;
         for (int j = 0; j < A->ig[num_of_knots]; j++)
            if (A->jg[j] == cond.knots_num[i])
               A->u[j] = 0.;
         //b[cond.knots_num[i]] = 1e10 * ug(mesh->knots[cond.knots_num[i]]);
         b[cond.knots_num[i]] = ug(mesh->knots[cond.knots_num[i]]);
      }
   }
#endif
}

void FEM::AddSecondBounds()
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
         h2 = pow(r2, 2.) - pow(r1, 2.);
         real h2_4 = h*h/4.,
              hr_3 = h*r1/3.;
         localMb[0][0] = hr_3      + h2_4 / 3.;  localMb[0][1] = hr_3 / 2. + h2_4 / 3.;
         localMb[1][0] = hr_3 / 2. + h2_4 / 3.;  localMb[1][1] = hr_3      + h2_4;
      }
      else
      {
         h = z2 - z1;
         h2 = pow(z2, 2.) - pow(z1, 2.);
         real rh_3 = r1 * h / 3.;
         localMb[0][0] = rh_3;       localMb[0][1] = rh_3 / 2.;
         localMb[1][0] = rh_3 / 2.;  localMb[1][1] = rh_3;
      }
      //real ratio = len;
      

      for (int i = 0; i < 2; i++)
         b[bound.knots_num[i]] += bound.value1 * (localMb[i][0] + localMb[i][1]) / (h2 * 3.1415926535897);
   }
}

void FEM::AddThirdBounds()
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
         h2 = pow(r2, 2.) - pow(r1, 2.);
         real h2_4 = h*h/4.,
              hr_3 = h*r1/3.;
         localMb[0][0] = hr_3      + h2_4 / 3.;  localMb[0][1] = hr_3 / 2. + h2_4 / 3.;
         localMb[1][0] = hr_3 / 2. + h2_4 / 3.;  localMb[1][1] = hr_3      + h2_4;
      }
      else
      {
         h = z2 - z1;
         h2 = pow(z2, 2.) - pow(z1, 2.);
         real rh_3 = r1 * h / 3.;
         localMb[0][0] = rh_3;       localMb[0][1] = rh_3 / 2.;
         localMb[1][0] = rh_3 / 2.;  localMb[1][1] = rh_3;
      }
      //real ratio = len;

      for (int i = 0; i < 2; i++)
      {
         b[bound.knots_num[i]] += bound.value1 * bound.value2 * (localMb[i][0] + localMb[i][1]) / (h2 * 3.1415926535897);
         for (int j = 0; j < 2; j++)
            AddElement(A, bound.knots_num[i], bound.knots_num[j], localMb[i][j]);
      }
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
      CreateExtraG(elem);
      Createb(elem);
      AddToGlobalMatricies(elem);
   }

   AssembleMatricies(isTimed);
   AddSecondBounds();
   AddThirdBounds();
   AddFirstBounds();
}

void FEM::AssembleMatricies(bool isTimed)
{
   real mulM = isTimed ? 1. / th : 1.;
   real mulGx = isTimed ? 0.0 : 1.;
   for (int i = 0; i < A->l.size(); i++)
   {
      A->l[i] = mulM * M->l[i] + G->l[i] + mulGx * Gx->l[i];
      A->u[i] = mulM * M->u[i] + G->u[i] + mulGx * Gx->u[i];
   }
   for (int i = 0; i < A->dim; i++)
      A->di[i] = mulM * M->di[i] + G->di[i] + mulGx * Gx->di[i];

   if (isTimed)
   {
      std::vector<real> Mq;
      Mq.resize(b.size());
      MatxVec(Mq, M, q1);
      for (int i = 0; i < Mq.size(); i++)
         d[i] = b[i] + Mq[i];
   }
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
   real ht2_12 = h * t * t / 12., 
        htr_9 = h * t * r1 / 9;

   localM[0][0] = localM[2][2] = ht2_12 / 3. + htr_9;
   localM[1][1] = localM[3][3] = ht2_12      + htr_9;
   localM[3][0] = localM[0][3] = 
   localM[1][2] = localM[2][1] = ht2_12 / 6. + htr_9 / 4.;
   localM[0][1] = localM[1][0] = 
   localM[2][3] = localM[3][2] = ht2_12 / 3. + htr_9 / 2.;
   localM[0][2] = localM[2][0] = ht2_12 / 6. + htr_9 / 2.;
   localM[1][3] = localM[3][1] = ht2_12 / 2. + htr_9 / 2.;

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         //localM[i][j] = Integrate(Mij, i, j, elem->knots_num);
         localM[i][j] *= elem.gam;
               
}

void FEM::CreateG(element& elem)
{
   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real h = r2 - r1, t = z2 - z1;
   real h2_4t = h * h / t / 4.,
        hr_3t = h * r1 / t / 3.,
        rt_3h = r1 * t / h / 3.;

   localG[0][0] = localG[2][2] =  h2_4t / 3. + hr_3t      + t / 6.  + rt_3h;
   localG[1][1] = localG[3][3] =  h2_4t      + hr_3t      + t / 6.  + rt_3h;
   localG[3][0] = localG[2][1] = 
   localG[1][2] = localG[0][3] = -h2_4t / 3. - hr_3t / 2. - t / 12. - rt_3h / 2.;
   localG[0][1] = localG[1][0] = 
   localG[3][2] = localG[2][3] =  h2_4t / 3. + hr_3t / 2. - t / 6.  - rt_3h;
   localG[1][3] = localG[3][1] = -h2_4t      - hr_3t      + t / 12. + rt_3h / 2.;
   localG[0][2] = localG[2][0] = -h2_4t / 3. - hr_3t      + t / 12. + rt_3h / 2.;

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localG[i][j] *= elem.lam;//* Integrate(Gij, i, j, elem->knots_num);

}

void FEM::CreateExtraG(element& elem)
{

   knot v1 = mesh->GetVelocity(mesh->knots[elem.knots_num[0]], 1);
   knot v2 = mesh->GetVelocity(mesh->knots[elem.knots_num[1]], 1);
   knot v3 = mesh->GetVelocity(mesh->knots[elem.knots_num[2]], 1);
   knot v4 = mesh->GetVelocity(mesh->knots[elem.knots_num[3]], 1);

   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real h = r2 - r1, t = z2 - z1;
   real ht2_9 = h * t * t / 9.,
        t2p_6 = t * t * r1 / 6.,
        t3_8 = t * t * t / 8.,
        h2p_6 = h * h * r1 / 6.;

   localGx[0][0] = v1.x * (-ht2_9 / 2. - t2p_6     ) + v1.y * (-t3_8 / 3. - h2p_6     );
   localGx[1][0] = v2.x * (-ht2_9      - t2p_6     ) + v1.y * (-t3_8 / 3. - h2p_6 / 2.);
   localGx[2][0] = v3.x * (-ht2_9 / 4. - t2p_6 / 2.) + v1.y * (-t3_8 / 3. - h2p_6     );
   localGx[3][0] = v4.x * (-ht2_9 / 2. - t2p_6 / 2.) + v1.y * (-t3_8 / 3. - h2p_6 / 2.);
                           
   localGx[0][1] = v1.x * (-ht2_9 / 2. - t2p_6     ) + v2.y * (t3_8 / 3. + h2p_6 / 2.);
   localGx[1][1] = v2.x * (-ht2_9      - t2p_6     ) + v2.y * (t3_8      + h2p_6     );
   localGx[2][1] = v3.x * (-ht2_9 / 4. - t2p_6 / 2.) + v2.y * (t3_8 / 3. + h2p_6 / 2.);
   localGx[3][1] = v4.x * (-ht2_9 / 2. - t2p_6 / 2.) + v2.y * (t3_8      + h2p_6     );
                           
   localGx[0][2] = v1.x * ( ht2_9 / 4. + t2p_6 / 2.) + v3.y * (-t3_8 / 3. - h2p_6);
   localGx[1][2] = v2.x * ( ht2_9 / 2. + t2p_6 / 2.) + v3.y * (-t3_8 / 3. - h2p_6 / 2.);
   localGx[2][2] = v3.x * ( ht2_9 / 2. + t2p_6     ) + v3.y * (-t3_8 / 3. - h2p_6);
   localGx[3][2] = v4.x * ( ht2_9 / 1. + t2p_6     ) + v3.y * (-t3_8 / 3. - h2p_6 / 2.);
                           
   localGx[0][3] = v1.x * ( ht2_9 / 4. + t2p_6 / 2.) + v4.y * (t3_8 / 3. + h2p_6 / 2.);
   localGx[1][3] = v2.x * ( ht2_9 / 2. + t2p_6 / 2.) + v4.y * (t3_8      + h2p_6     );
   localGx[2][3] = v3.x * ( ht2_9 / 2. + t2p_6     ) + v4.y * (t3_8 / 3. + h2p_6 / 2.);
   localGx[3][3] = v4.x * ( ht2_9 / 1. + t2p_6     ) + v4.y * (t3_8      + h2p_6     );

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localGx[i][j] *= elem.gam;//* Integrate(Gij, i, j, elem->knots_num);
}

void FEM::Createb(element& elem)
{
   real localb[4]{};
   real f_[4]{};
   for (int i = 0; i < 4; i++)
      f_[i] = f(mesh->knots[elem.knots_num[i]]);

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localb[i] += localM[i][j] * f_[j];

   for (int i = 0; i < 4; i++)
      b[elem.knots_num[i]] += localb[i] / elem.gam;
}
