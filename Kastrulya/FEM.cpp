#include "FEM.h"
#include <string>
FEM::FEM(Mesh* _mesh)
{
   //std::ifstream fknots("Knots.txt");
   //std::ifstream felems("elements.txt");
   //std::ifstream fbounds1("FirstBounds.txt");
   //std::ifstream fbounds2("SecondBounds.txt");
   //std::ifstream fbounds3("ThirdBounds.txt");
   std::ifstream fparams("EnvParams.txt");
   fparams >> un;
   fparams.close();

   mesh->MakeMesh();

   A = MakeSparseFormat(4, mesh->knots.size(), mesh);
   M = MakeSparseFormat(4, mesh->knots.size(), mesh);
   G = MakeSparseFormat(4, mesh->knots.size(), mesh);

   q.resize(num_of_knots, 0.);
   b.resize(num_of_knots, 0.);

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
   CreateSLAE();
   SolveSLAE(A, q, b);

#ifdef DEBUG
   check_test();
#endif // DEBUG
}

void FEM::SolveParabolic()
{
   for  (real t = 0, int tn = 0; t < t_last; t += th, tn++)
   {
      if (tn > 1)
      {
         CreateSLAE();
         SolveSLAE(A, q, b);
      }
      std::string str = "./Results/Result_layer_" + std::to_string(tn) + ".txt";
      std::ofstream out(str);
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
   out << "z";
   out.width(15);
   out << "q";
   out.width(15);
   out << "\n";
   out << std::setprecision(14);

   for (int i = 0; i < num_of_knots; i++)
   {
      out << std::defaultfloat;
      out << mesh->knots[i].x;
      out.width(15);
      out << mesh->knots[i].y;
      out.width(15);
      out << mesh->knots[i].z;
      out.width(15);
      out << std::scientific;
      out << q[i];
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
      real len = 0;

      real ratio = len;
      localMb[0][0] = len / 3.;  localMb[0][1] = len / 6.;
      localMb[1][0] = len / 6.;  localMb[1][1] = len / 3.;

      for (int i = 0; i < 2; i++)
         b[bound.knots_num[i]] += bound.value1 * (localMb[i][0] + localMb[i][1]) / len;
   }
}

void FEM::CreateSLAE()
{
   element elem;
   for (int i = 0; i < num_of_FE; i++)
   {
      elem = mesh->elems[i];
      CreateG(elem);
      CreateM(elem);
      //AddToA(elem);
      AssembleMatricies(elem);
      //AddFirstBounds();
   }
   AddSecondBounds();
   AddThirdBounds();
}

void FEM::AssembleMatricies(element& elem)
{
   for (int i = 0; i < A->l.size(); i++)
   {
      A->l[i] = 1. / th * elem.gam * M->l[i] + G->l[i];
      A->u[i] = 1. / th * elem.gam * M->u[i] + G->u[i];
   }
   for (int i = 0; i < A->dim; i++)
      A->di[i] = 1. / th * elem.gam * M->di[i] + G->di[i];
}

void FEM::AddToA(element& elem)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
      {
         localA[i][j] = localG[i][j] + localM[i][j];
         AddElement(A, elem.knots_num[i], elem.knots_num[j], localA[i][j]);
      }
//AddLocal(A, elem->knots_num, localA, 1);
}

void FEM::CreateM(element& elem)
{
   real x1 = mesh->knots[elem.knots_num[0]].x;
   real x2 = mesh->knots[elem.knots_num[1]].x;
   real y1 = mesh->knots[elem.knots_num[0]].y;
   real y2 = mesh->knots[elem.knots_num[2]].y;
   real hx = x2 - x1, hy = y2 - y1, c = hx * hy / 36.;

   localM[0][0] = localM[1][1] = localM[2][2] = localM[3][3] = 4. * c;
   localM[0][3] = localM[1][2] = localM[2][1] = localM[3][0] = c;
   localM[0][1] = localM[0][2] = localM[1][0] = localM[1][3] = 
   localM[2][0] = localM[2][3] = localM[3][1] = localM[3][2] = 2. * c;

   //for (int i = 0; i < 8; i++)
   //   for (int j = 0; j < 8; j++)
   //      localM[i][j] = Integrate(Mij, i, j, elem->knots_num);
}

void FEM::CreateG(element& elem)
{
   real r1 = mesh->knots[elem.knots_num[0]].x;
   real r2 = mesh->knots[elem.knots_num[1]].x;
   real z1 = mesh->knots[elem.knots_num[0]].y;
   real z2 = mesh->knots[elem.knots_num[2]].y;
   real h = r2 - r1, t = z2 - z1;
   

   localG[0][0] = h * h / 12. / t + h * r1 / 3. / t + t / 6. + r1 * t / 3 / h; 
   localG[1][1] = 
   localG[2][2] = localG[3][3] = 2. * ( cyx + cxy );
   localG[3][0] = localG[2][1] = localG[1][2] = localG[0][3] = - cyx - cxy;
   localG[0][1] = localG[1][0] = localG[3][2] = localG[2][3] = -2. * cyx + cxy;
   localG[1][3] = localG[3][1] = localG[0][2] = localG[2][0] = cyx - 2. * cxy;

   //for (int i = 0; i < 4; i++)
   //   for (int j = 0; j < 4; j++)
   //      localG[i][j] = elem.lam * Integrate(Gij, i, j, elem->knots_num);

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
      b[elem.knots_num[i]] += localb[i];
}
