#include "./roots/roots_solver.h"
#include <stdio.h>

double f(double x)
{
    return x * x - 2.0;
}

void logger(const char *log)
{
    printf("%s\n", log);
}

int main()
{
    RootSolver rs = rootsolver_create(f, NULL, logger);

    double sol;
    RootSolverStatus status = rootsolver_bisection(&rs, &sol, 0.0, 2.0, 5e-9, 100);

    return 0;
}