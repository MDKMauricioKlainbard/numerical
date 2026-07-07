#include "ivp_solver.h"
#include <stdio.h>

IVPSolver ivpsolver_create(OdeFunc f, int n, LogFunc log) {
    IVPSolver ivp;
    ivp.f = f;
    ivp.n = n;
    ivp.log = log;
    return ivp;
}

void ivpsolver_step_euler(const IVPSolver *ivp, double t, double *y, double dt) {
    
}