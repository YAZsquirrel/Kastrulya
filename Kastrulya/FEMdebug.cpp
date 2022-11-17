#include "FEM.h"

#ifdef DEBUG
real FEM::f(knot& k, element& e)
{
	switch (e.n_mat)
	{
		case 1: return 2. + k.y * k.y;
		case 2: return -6. * k.y + k.y * k.y * k.y;
		case 3: return k.y;
		case 4: return 0;
		case 5: return k.x * k.x * k.x * k.x + k.y * k.y + k.x * k.y + 1 - 12 * k.x * k.x - 2;

	default:
		return 0;
	}
}

real FEM::bound1func(knot& k, int n_mat) 
{
	switch (n_mat)
	{
	case 1: return 2. + k.y * k.y;
	case 2: return -6. * k.y + k.y * k.y * k.y;
	case 3: return k.y;
	case 4: return 0;
	case 5: return k.x * k.x * k.x * k.x + k.y * k.y + k.x * k.y + 1 - 12 * k.x * k.x - 2;

	default:
		return 0;
	}
}

real FEM::bound2func(knot& k, int n_mat)
{
	switch (n_mat)
	{
	case 1: return 2. + k.y * k.y;
	case 2: return -6. * k.y + k.y * k.y * k.y;
	case 3: return k.y;
	case 4: return 0;
	case 5: return k.x * k.x * k.x * k.x + k.y * k.y + k.x * k.y + 1 - 12 * k.x * k.x - 2;

	default:
		return 0;
	}
}

real FEM::bound3func(knot& k, int n_mat)
{
	switch (n_mat)
	{
	case 1: return 2. + k.y * k.y;
	case 2: return -6. * k.y + k.y * k.y * k.y;
	case 3: return k.y;
	case 4: return 0;
	case 5: return k.x * k.x * k.x * k.x + k.y * k.y + k.x * k.y + 1 - 12 * k.x * k.x - 2;

	default:
		return 0;
	}
}

real FEM::bound3funcbeta(knot& k, int n_mat)
{
	switch (n_mat)
	{
	case 1: return 2. + k.y * k.y;
	case 2: return -6. * k.y + k.y * k.y * k.y;
	case 3: return k.y;
	case 4: return 0;
	case 5: return k.x * k.x * k.x * k.x + k.y * k.y + k.x * k.y + 1 - 12 * k.x * k.x - 2;

	default:
		return 0;
	}
}
#endif // DEBUG
