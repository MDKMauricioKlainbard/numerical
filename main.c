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
    RootSolver rs = rootsolver_create(f, logger); // <-- f, SIN paréntesis ni &

    double sol;
    RootSolverStatus status = rootsolver_bisection(&rs, &sol, 0.0, 2.0, 1e-9, 100);

    if (status == ROOTSOLVER_OK)
    {
        printf("Solucion encontrada: %.10f\n", sol);
    }
    else
    {
        printf("Error: codigo %d\n", status);
    }

    return 0;
}