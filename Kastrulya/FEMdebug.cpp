#include "FEM.h"
#include <math.h>

#ifdef DEBUG2
real FEM::f(knot& k, element& e, real t)
{
   switch (e.n_test)
   {
      case 0:
         if (abs(k.x) > 1e-12)
            return -1 / k.x;		// check G
         else return 0;

      case 1: return k.x;	// check M
      case 2: return k.x;	// check G + M
      case 3: return k.y;
      case 4: return 51;
      case 5:
         if (abs(k.x) > 1e-12)
            return k.x * k.y - k.y / k.x;  // G + M
         else return k.x * k.y;
      case 6: return -4;            // G 
      case 7: return k.x * k.x;     // M
      case 8: return k.x * k.x - 4; // G + M
      case 9: 
         if (abs(k.x) > 1e-12)
            return 2 * sin(k.x) - cos(k.x) / k.x;  // G + M
         else return sin(k.x);
      case 10: return sin(k.x);  // M

      // t ->
      case 11:
         if (abs(k.x) > 1e-12)
            return k.x - t / k.x;  // G + M
         else return k.x;
      case 12:
         if (abs(k.x) > 1e-12)
            return - 1. / k.x;  // G + M
         else return 0;
      case 13:
         return 1.; // 1 ????     // G + M
      case 14:
         return 1.; // 1 ????     // G + M
      case 15:
         return 0.; // 1 ????     // G + M
   default:
      return 0;
   }
}

real FEM::bound1func(knot& k, real t, int n_test)
{
   switch (n_test)
   {
      case 0: return k.x;  // r
      case 1: return k.x;  // r
      case 2: return k.x;	// r
      case 3: return k.y;  // z
      case 4: return 51;   // 51
      case 5: return k.x * k.y; // xy 
      case 6: return k.x * k.x; // -divgrad(x**2)
      case 7: return k.x * k.x; //  x**2
      case 8: return k.x * k.x; // -divgrad(x**2) + x**2
      case 9: return sin(k.x);
      case 10: return sin(k.x);

      // t ->
      case 11: return k.x * t;
      case 12: return k.x;
      case 13: return t;
      case 14: return k.y + t;      // +
      case 15: return k.y;      // +

   default:
      return 0;
   }
}

real FEM::bound2func(knot& k, real time, int n_test)
{
   switch (n_test)
   {
   case 5: return k.y;

   default:
      return 0;
   }
}

real FEM::bound3func(knot& k, real time, int n_test)
{
   switch (n_test)
   {
   case 1: return 2. + k.y * k.y;

   default:
      return 0;
   }
}

real FEM::bound3funcbeta(knot& k, real time, int n_test)
{
   switch (n_test)
   {
   case 1: return 2. + k.y * k.y;

   default:
      return 0;
   }
}
#endif // DEBUG
