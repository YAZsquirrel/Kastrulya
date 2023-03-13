#include "FEM.h"
#include <mkl_pardiso.h>


int main()
//int _tmain(int argc, char* argv[])
{
	FEM *fem = new FEM();
	fem->SolveParabolic();
	return 0;
}