#include "ivp_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

static void log_msg(const IVPSolver *ivp, const char *fmt, ...)
{
    if (!ivp->log)
        return;

    char buf[128];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);

    ivp->log(buf);
}

IVPSolverStatus ivpsolver_create(IVPSolver *ivp, OdeFunc f, int n, LogFunc log)
{
    ivp->f = f;
    ivp->n = n;
    ivp->log = log;
    ivp->dydt_buffer = malloc(n * sizeof(double));

    if (ivp->dydt_buffer == NULL)
    {
        log_msg(ivp, "[IVP SOLVER] - Memory allocation failed during creation.");
        return IVPSOLVER_ALLOC_FAILED;
    }

    return IVPSOLVER_OK;
}

void ivpsolver_destroy(IVPSolver *ivp)
{
    free(ivp->dydt_buffer);
    ivp->dydt_buffer = NULL;
}

void ivpsolver_step_euler(const IVPSolver *ivp, double t, double *y, double dt)
{

    ivp->f(t, y, ivp->dydt_buffer, ivp->n);

    for (int i = 0; i < ivp->n; i++)
    {
        y[i] = y[i] + ivp->dydt_buffer[i] * dt;
    }
}