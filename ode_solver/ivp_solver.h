#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H

typedef void (*OdeFunc)(double t, const double *y, double *dydt, int n);
typedef void (*LogFunc)(const char *msg);

typedef struct
{
    OdeFunc f;
    int n; // Dimension
    LogFunc log;
} IVPSolver;

IVPSolver ivpsolver_create(OdeFunc f, int n, LogFunc log);

void ivpsolver_step_euler(const IVPSolver *ivp, double t, double *y, double dt);

#endif