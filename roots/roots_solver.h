#ifndef ROOTS_SOLVER_H
#define ROOTS_SOLVER_H

typedef double (*FuncPtr)(double);
typedef void (*LogFunc)(const char* msg);

typedef struct {
    FuncPtr f;
    LogFunc log;
} RootSolver;

typedef enum {
    ROOTSOLVER_OK = 0,
    ROOTSOLVER_INVALID_INTERVAL = 1,
    ROOTSOLVER_NO_CONVERGENCE = 2
} RootSolverStatus;

RootSolver rootsolver_create(FuncPtr f, LogFunc log);

// Closed Methods
RootSolverStatus rootsolver_bisection(const RootSolver* rs, double *sol, double a, double b, double tol, int max_iter);


#endif