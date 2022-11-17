#include "FEM.h"


int main()
{
	FEM *fem = new FEM();
	fem->SolveParabolic();
	return 0;
}