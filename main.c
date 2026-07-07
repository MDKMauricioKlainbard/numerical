#include "./roots/roots_solver.h"
#include <stdio.h>
#include <math.h>

// f(x) = x^2 - 2, raiz real: sqrt(2)
double f(double x)
{
    return x * x - 2.0;
}

// f'(x) = 2x, derivada analitica para Newton
double df(double x)
{
    return 2.0 * x;
}

// g(x) = (x + 2/x) / 2, reescritura de f(x)=0 como x = g(x)
// Este despeje viene del metodo de Newton aplicado a x^2 = 2
double g(double x)
{
    return (x + 2.0 / x) / 2.0;
}

// Logger simple: imprime cada mensaje interno del solver
void console_logger(const char *msg)
{
    printf("%s\n", msg);
}

static void print_result(const char *method_name, RootSolverStatus status, double sol)
{
    printf("\n--- %s ---\n", method_name);
    if (status == ROOTSOLVER_OK)
    {
        printf("Resultado: %.10f\n", sol);
        printf("Error vs sqrt(2): %.2e\n", fabs(sol - sqrt(2.0)));
    }
    else
    {
        printf("Fallo con codigo de estado: %d\n", status);
    }
}

int main()
{
    double sol;
    RootSolverStatus status;

    double tol = 1e-9;
    int max_iter = 100;

    // --- Bisection: necesita f y un intervalo que encierre la raiz ---
    RootSolver rs_bisection = rootsolver_create(f, NULL, console_logger);
    status = rootsolver_bisection(&rs_bisection, &sol, 0.0, 2.0, tol, max_iter);
    print_result("BISECTION", status, sol);

    // --- Regula Falsi: mismo intervalo que bisection ---
    RootSolver rs_regulafalsi = rootsolver_create(f, NULL, console_logger);
    status = rootsolver_regulafalsi(&rs_regulafalsi, &sol, 0.0, 2.0, tol, max_iter);
    print_result("REGULA FALSI", status, sol);

    // --- Fixed Point: usa g(x), no f(x) directamente ---
    RootSolver rs_fixedpoint = rootsolver_create(g, NULL, console_logger);
    status = rootsolver_fixedpoint(&rs_fixedpoint, &sol, 1.0, tol, max_iter);
    print_result("FIXED POINT", status, sol);

    // --- Newton-Raphson con derivada analitica explicita ---
    RootSolver rs_newton = rootsolver_create(f, df, console_logger);
    status = rootsolver_newtonraphson(&rs_newton, &sol, 1.0, tol, max_iter);
    print_result("NEWTON (con df explicita)", status, sol);

    // --- Newton-Raphson SIN derivada explicita: usa diferencias finitas ---
    RootSolver rs_newton_numeric = rootsolver_create(f, NULL, console_logger);
    status = rootsolver_newtonraphson(&rs_newton_numeric, &sol, 1.0, tol, max_iter);
    print_result("NEWTON (diferencias finitas)", status, sol);

    status = rootsolver_bisection(&rs_bisection, &sol, 3.0, 5.0, tol, max_iter);
    print_result("BISECTION (intervalo invalido)", status, sol);

    return 0;
}