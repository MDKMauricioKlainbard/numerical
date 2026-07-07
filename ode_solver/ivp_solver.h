#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H

typedef void (*OdeFunc)(double t, const double *y, double *dydt, int n);
typedef void (*LogFunc)(const char *msg);

typedef struct
{
    OdeFunc f;
    int n; // Dimension
    LogFunc log;
    double *dydt_buffer;
} IVPSolver;

typedef enum
{
    IVPSOLVER_OK = 0,
    IVPSOLVER_ALLOC_FAILED = 1
} IVPSolverStatus;

IVPSolverStatus ivpsolver_create(IVPSolver *ivp, OdeFunc f, int n, LogFunc log);

void ivpsolver_step_euler(const IVPSolver *ivp, double t, double *y, double dt);
void ivpsolver_step_rk4(const IVPSolver *ivp, double t, double *y, double dt);

void ivpsolver_destroy(IVPSolver *ivp);

#endif