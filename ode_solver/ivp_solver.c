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
    ivp->k1 = malloc(n * sizeof(double));
    ivp->k2 = malloc(n * sizeof(double));
    ivp->k3 = malloc(n * sizeof(double));
    ivp->k4 = malloc(n * sizeof(double));
    ivp->y_temp = malloc(n * sizeof(double));

    if (!ivp->k1 || !ivp->k2 || !ivp->k3 || !ivp->k4 || !ivp->y_temp)
    {
        log_msg(ivp, "[IVP SOLVER] - Memory allocation failed during creation.");
        return IVPSOLVER_ALLOC_FAILED;
    }

    return IVPSOLVER_OK;
}

void ivpsolver_destroy(IVPSolver *ivp)
{
    free(ivp->k1);
    free(ivp->k2);
    free(ivp->k3);
    free(ivp->k4);
    free(ivp->y_temp);
    ivp->k1 = ivp->k2 = ivp->k3 = ivp->k4 = ivp->y_temp = NULL;
}

void ivpsolver_step_euler(const IVPSolver *ivp, double t, double *y, double dt)
{
    int n = ivp->n;

    ivp->f(t, y, ivp->k1, n);

    for (int i = 0; i < n; i++)
    {
        y[i] = y[i] + ivp->k1[i] * dt;
    }
}

void ivpsolver_step_rk4(const IVPSolver *ivp, double t, double *y, double dt)
{
    int n = ivp->n;

    ivp->f(t, y, ivp->k1, n);

    for (int i = 0; i < n; i++)
        ivp->y_temp[i] = y[i] + (dt / 2.0) * ivp->k1[i];
    ivp->f(t + dt / 2.0, ivp->y_temp, ivp->k2, n);

    for (int i = 0; i < n; i++)
        ivp->y_temp[i] = y[i] + (dt / 2.0) * ivp->k2[i];
    ivp->f(t + dt / 2.0, ivp->y_temp, ivp->k3, n);

    for (int i = 0; i < n; i++)
        ivp->y_temp[i] = y[i] + dt * ivp->k3[i];
    ivp->f(t + dt, ivp->y_temp, ivp->k4, n);

    for (int i = 0; i < n; i++)
        y[i] = y[i] + (dt / 6.0) * (ivp->k1[i] + 2.0 * ivp->k2[i] + 2.0 * ivp->k3[i] + ivp->k4[i]);
}