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
         SolveSLAE_LOS(A, q1, d);

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
}

void FEM::AddFirstBounds(real time)
{
   for (auto& cond : mesh->bounds1)
   {
      for (int i = 0; i < 2; i++)
      {
         A->di[cond.knots_num[i]] = 1.;
         for (int j = A->ig[cond.knots_num[i]]; j < A->ig[cond.knots_num[i] + 1]; j++)
            A->l[j] = 0.;
         for (int ii = 0; ii < A->dim; ii++)                // идем по столбцам
            for (int j = A->ig[ii]; j < A->ig[ii + 1]; j++)   // идем элементам в столбце
               if (A->jg[j] == cond.knots_num[i])          // в нужной строке элемент?
                  A->u[j] = 0.;

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

   localC[0][0] = 
   localC[1][1] = h2t_72 * 2. + htr_36 * 4.;       // +
   localC[2][2] = 
   localC[3][3] = h2t_72 * 6. + htr_36 * 4.;       // +
   localC[0][3] = localC[3][0] =                            // +
   localC[1][2] = localC[2][1] = h2t_72      + htr_36;  // +
   localC[0][2] = localC[2][0] =                            //
   localC[1][3] = localC[3][1] = h2t_72 * 2. + htr_36 * 2.;  //
   localC[0][1] = localC[1][0] = h2t_72      + htr_36 * 2.;  // +
   localC[2][3] = localC[3][2] = h2t_72 * 3. + htr_36 * 2.;  // +

   // 2 1 2 1
   // 1 2 1 2
   // 2 1 6 3
   // 1 2 3 6

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         //localM[i][j] = Integrate(Mij, i, j, elem->knots_num);
         localM[i][j] = elem.gam * localC[i][j];
               
}

void FEM::CreateG(element& elem)
{
   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real h = r2 - r1, t = z2 - z1;
   //real h2_4t = h * h / t / 4.,
   //     hr_3t = h * r1 / t / 3.,
   //     rt_3h = r1 * t / h / 3.;
   //
   //localG[0][0] = localG[2][2] =  h2_4t / 3. + hr_3t      + t / 6.  + rt_3h;
   //localG[1][1] = localG[3][3] =  h2_4t      + hr_3t      + t / 6.  + rt_3h;
   //localG[3][0] = localG[2][1] = 
   //localG[1][2] = localG[0][3] = -h2_4t / 3. - hr_3t / 2. - t / 12. - rt_3h / 2.;
   //localG[0][1] = localG[1][0] = 
   //localG[3][2] = localG[2][3] =  h2_4t / 3. + hr_3t / 2. - t / 6.  - rt_3h;
   //localG[1][3] = localG[3][1] = -h2_4t      - hr_3t      + t / 12. + rt_3h / 2.;
   //localG[0][2] = localG[2][0] = -h2_4t / 3. - hr_3t      + t / 12. + rt_3h / 2.;
   real ht2r_12hk = t * (2 * r1 + h) / h / 12.,
        hr_6t = h * r1 / t / 6.,
        h2_12t = h * h / t / 12.;
   localG[0][0] = 
   localG[2][2] = 2. * ht2r_12hk + 2. * hr_6t +      h2_12t;
   localG[1][1] = 
   localG[3][3] = 2. * ht2r_12hk + 2. * hr_6t + 3. * h2_12t;
   localG[1][2] = localG[2][1] = 
   localG[0][3] = localG[3][0] = - ht2r_12hk - hr_6t - h2_12t;
   localG[0][1] = localG[1][0] =
   localG[2][3] = localG[3][2] =  -2 * ht2r_12hk + hr_6t + h2_12t;
   localG[1][3] = localG[3][1] =  ht2r_12hk - 2. * hr_6t - 3. * h2_12t;
   localG[0][2] = localG[2][0] =  ht2r_12hk - 2. * hr_6t - h2_12t;



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
   real ht_9 = h * t / 9.,
        tp_6 = t  * r1 / 6.,
        t2_8 = t * t / 8.,
        hp_6 = h * r1 / 6.;

   //localGx[0][0] = v1.x * (-ht_9 / 2. - tp_6     ) + v1.y * (-t2_8 / 3. - hp_6     );
   //localGx[1][0] = v2.x * (-ht_9      - tp_6     ) + v1.y * (-t2_8 / 3. - hp_6 / 2.);
   //localGx[2][0] = v3.x * (-ht_9 / 4. - tp_6 / 2.) + v1.y * (-t2_8 / 3. - hp_6     );
   //localGx[3][0] = v4.x * (-ht_9 / 2. - tp_6 / 2.) + v1.y * (-t2_8 / 3. - hp_6 / 2.);
   //                        
   //localGx[0][1] = v1.x * (-ht_9 / 2. - tp_6     ) + v2.y * (t2_8 / 3. + hp_6 / 2.);
   //localGx[1][1] = v2.x * (-ht_9      - tp_6     ) + v2.y * (t2_8      + hp_6     );
   //localGx[2][1] = v3.x * (-ht_9 / 4. - tp_6 / 2.) + v2.y * (t2_8 / 3. + hp_6 / 2.);
   //localGx[3][1] = v4.x * (-ht_9 / 2. - tp_6 / 2.) + v2.y * (t2_8      + hp_6     );
   //                        
   //localGx[0][2] = v1.x * ( ht_9 / 4. + tp_6 / 2.) + v3.y * (-t2_8 / 3. - hp_6);
   //localGx[1][2] = v2.x * ( ht_9 / 2. + tp_6 / 2.) + v3.y * (-t2_8 / 3. - hp_6 / 2.);
   //localGx[2][2] = v3.x * ( ht_9 / 2. + tp_6     ) + v3.y * (-t2_8 / 3. - hp_6);
   //localGx[3][2] = v4.x * ( ht_9 / 1. + tp_6     ) + v3.y * (-t2_8 / 3. - hp_6 / 2.);
   //                        
   //localGx[0][3] = v1.x * ( ht_9 / 4. + tp_6 / 2.) + v4.y * (t2_8 / 3. + hp_6 / 2.);
   //localGx[1][3] = v2.x * ( ht_9 / 2. + tp_6 / 2.) + v4.y * (t2_8      + hp_6     );
   //localGx[2][3] = v3.x * ( ht_9 / 2. + tp_6     ) + v4.y * (t2_8 / 3. + hp_6 / 2.);
   //localGx[3][3] = v4.x * ( ht_9 / 1. + tp_6     ) + v4.y * (t2_8      + hp_6     ); 

   localGx[0][0] = v1.x * (-ht_9 / 2. - tp_6     ) + v1.y * (-t2_8 / 3. - hp_6     );
   localGx[1][0] = v2.x * (-ht_9      - tp_6     ) + v1.y * (-t2_8 / 3. - hp_6 / 2.);
   localGx[2][0] = v3.x * (-ht_9 / 4. - tp_6 / 2.) + v1.y * (-t2_8 / 3. - hp_6     );
   localGx[3][0] = v4.x * (-ht_9 / 2. - tp_6 / 2.) + v1.y * (-t2_8 / 3. - hp_6 / 2.);
   
   localGx[0][1] = v1.x * ( ht_9 / 2. + tp_6     ) + v2.y * (-t2_8 / 3. - hp_6 / 2.);
   localGx[1][1] = v2.x * ( ht_9      + tp_6     ) + v2.y * (-t2_8      - hp_6     );
   localGx[2][1] = v3.x * ( ht_9 / 4. + tp_6 / 2.) + v2.y * (-t2_8 / 3. - hp_6 / 2.);
   localGx[3][1] = v4.x * ( ht_9 / 2. + tp_6 / 2.) + v2.y * (-t2_8      - hp_6     );
   
   localGx[0][2] = v1.x * (-ht_9 / 4. - tp_6 / 2.) + v3.y * ( t2_8 / 3. + hp_6     );
   localGx[1][2] = v2.x * (-ht_9 / 2. - tp_6 / 2.) + v3.y * ( t2_8 / 3. + hp_6 / 2.);
   localGx[2][2] = v3.x * (-ht_9 / 2. - tp_6     ) + v3.y * ( t2_8 / 3. + hp_6     );
   localGx[3][2] = v4.x * (-ht_9 / 1. - tp_6     ) + v3.y * ( t2_8 / 3. + hp_6 / 2.);
   
   localGx[0][3] = v1.x * ( ht_9 / 4. + tp_6 / 2.) + v4.y * ( t2_8 / 3. + hp_6 / 2.);
   localGx[1][3] = v2.x * ( ht_9 / 2. + tp_6 / 2.) + v4.y * ( t2_8      + hp_6     );
   localGx[2][3] = v3.x * ( ht_9 / 2. + tp_6     ) + v4.y * ( t2_8 / 3. + hp_6 / 2.);
   localGx[3][3] = v4.x * ( ht_9 / 1. + tp_6     ) + v4.y * ( t2_8      + hp_6     );

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         localGx[i][j] *= elem.gam;
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
