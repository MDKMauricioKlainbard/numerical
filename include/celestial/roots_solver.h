#ifndef ROOTS_SOLVER_H
#define ROOTS_SOLVER_H

typedef double (*FuncPtr)(double);
typedef void (*LogFunc)(const char *msg);

typedef struct
{
    FuncPtr f;
    FuncPtr df;
    LogFunc log;
} RootSolver;

typedef enum
{
    ROOTSOLVER_OK = 0,
    ROOTSOLVER_INVALID_INTERVAL = 1,
    ROOTSOLVER_NO_CONVERGENCE = 2
} RootSolverStatus;

RootSolver rootsolver_create(FuncPtr f, FuncPtr df, LogFunc log);

// Closed Methods
RootSolverStatus rootsolver_bisection(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter);
RootSolverStatus rootsolver_regulafalsi(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter);

// Open Methods
RootSolverStatus rootsolver_fixedpoint(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter);
RootSolverStatus rootsolver_newtonraphson(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter);
#endif