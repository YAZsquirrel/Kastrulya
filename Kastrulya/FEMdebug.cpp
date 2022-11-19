#include "FEM.h"

#ifdef DEBUG2
real FEM::f(knot& k, element& e, real time)
{
	switch (e.n_test)
	{
		case 0: return 0;		// check G
		case 1: return k.x;	// check M
		case 2: return k.x;	// check G + M
		case 3: return k.y;

	default:
		return 0;
	}
}

real FEM::bound1func(knot& k, real time, int n_test)
{
	switch (n_test)
	{
	case 0: return k.x;
	case 1: return k.x;
	case 2: return k.x;
		case 3: return k.y;

	default:
		return 0;
	}
}

real FEM::bound2func(knot& k, real time, int n_test)
{
	switch (n_test)
	{
	case 0: return 2. + k.y * k.y;

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
