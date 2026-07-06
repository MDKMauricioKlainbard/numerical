#include "roots_solver.h"
#include <math.h>
#include <stdio.h>
#include <stdarg.h>

#define SAFETY_EPS 1e-12

static void log_msg(const RootSolver *rs, const char *fmt, ...)
{
    if (!rs->log)
        return;

    char buf[128];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);

    rs->log(buf);
}

RootSolver rootsolver_create(FuncPtr f, LogFunc log)
{
    RootSolver rs;
    rs.f = f;
    rs.log = log;
    return rs;
}

RootSolverStatus rootsolver_bisection(const RootSolver *rs, double *sol, double a, double b, double tol, int max_iter)
{
    if (a >= b)
    {
        log_msg(rs, "[BISECTION] - Invalid interval: a must be less than b");
        return ROOTSOLVER_INVALID_INTERVAL;
    }

    double fa = rs->f(a);
    double fb = rs->f(b);

    if ((fa * fb) > 0.0)
    {
        log_msg(rs, "[BISECTION] - No solution in interval [%.10f, %.10f]", a, b);
        return ROOTSOLVER_INVALID_INTERVAL;
    }

    int i = 1;
    double approx;
    double approx_old = a;
    double fapprox;

    while (i <= max_iter)
    {
        approx = (a + b) / 2.0;
        fapprox = rs->f(approx);
        double err_rel = fabs(approx - approx_old) / (fabs(approx) + SAFETY_EPS);

        if (err_rel < tol)
        {
            log_msg(rs, "[BISECTION] - Solution found: %.10f", approx);
            *sol = approx;
            return ROOTSOLVER_OK;
        }

        approx_old = approx;
        i += 1;

        if (fa * fapprox > 0.0)
        {
            a = approx;
            fa = fapprox;
        }
        else
        {
            b = approx;
            fb = fapprox;
        }
    }

    log_msg(rs, "[BISECTION] - Method failed after %d iterations. The procedure was unsuccessful.", max_iter);
    return ROOTSOLVER_NO_CONVERGENCE;
}