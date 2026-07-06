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

static double compute_derivative(const RootSolver *rs, double x)
{
    if (rs->df != NULL)
    {
        return rs->df(x);
    }
    else
    {
        double h = 1e-6;
        return (rs->f(x + h) - rs->f(x)) / h;
    }
}

RootSolver rootsolver_create(FuncPtr f, FuncPtr df, LogFunc log)
{
    RootSolver rs;
    rs.f = f;
    rs.df = df;
    rs.log = log;
    return rs;
}

RootSolverStatus rootsolver_bisection(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter)
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
    double err_rel;

    while (i <= max_iter)
    {
        approx = (a + b) / 2.0;
        fapprox = rs->f(approx);
        err_rel = fabs(approx - approx_old) / (fabs(approx) + SAFETY_EPS) * 100.0;

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

RootSolverStatus rootsolver_regulafalsi(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter)
{
    if (a >= b)
    {
        log_msg(rs, "[REGULA FALSI] - Invalid interval: a must be less than b");
        return ROOTSOLVER_INVALID_INTERVAL;
    }

    double fa = rs->f(a);
    double fb = rs->f(b);

    if ((fa * fb) > 0.0)
    {
        log_msg(rs, "[REGULA FALSI] - No solution in interval [%.10f, %.10f]", a, b);
        return ROOTSOLVER_INVALID_INTERVAL;
    }

    int i = 1;
    double approx;
    double approx_old = a;
    double fapprox;
    double err_rel;

    int left_repeat = 0;
    int right_repeat = 0;

    while (i <= max_iter)
    {
        approx = b - (fb * (a - b)) / (fa - fb);
        fapprox = rs->f(approx);
        err_rel = fabs(approx - approx_old) / (fabs(approx) + SAFETY_EPS) * 100.0;

        if (err_rel < tol)
        {
            log_msg(rs, "[REGULA FALSI] - Solution found: %.10f", approx);
            *sol = approx;
            return ROOTSOLVER_OK;
        }

        approx_old = approx;
        i += 1;

        if (fa * fapprox > 0.0)
        {
            a = approx;
            fa = fapprox;
            right_repeat += 1;
            left_repeat = 0;

            if (right_repeat >= 2)
            {
                fb /= 2.0;
            }
        }
        else
        {
            b = approx;
            fb = fapprox;
            left_repeat += 1;
            right_repeat = 0;

            if (left_repeat >= 2)
            {
                fa /= 2.0;
            }
        }
    }

    log_msg(rs, "[REGULA FALSI] - Method failed after %d iterations. The procedure was unsuccessful.", max_iter);
    return ROOTSOLVER_NO_CONVERGENCE;
}

RootSolverStatus rootsolver_fixedpoint(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter)
{
    double approx_old = approx;
    int i = 1;
    double err_rel;

    while (i <= max_iter)
    {
        approx = rs->f(approx_old);
        err_rel = fabs(approx - approx_old) / (fabs(approx) + SAFETY_EPS) * 100.0;

        if (err_rel < tol)
        {
            log_msg(rs, "[FIXED POINT] - Solution found: %.10f", approx);
            *sol = approx;
            return ROOTSOLVER_OK;
        }

        approx_old = approx;
        i += 1;
    }

    log_msg(rs, "[FIXED POINT] - Method failed after %d iterations. The procedure was unsuccessful.", max_iter);
    return ROOTSOLVER_NO_CONVERGENCE;
}

RootSolverStatus rootsolver_newtonraphson(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter)
{
    int i = 1;
    double approx_old = approx;
    double err_rel;
    double fx;
    double dfx;

    while (i <= max_iter)
    {
        fx = rs->f(approx_old);
        dfx = compute_derivative(rs, approx_old);

        if (fabs(dfx) < SAFETY_EPS)
        {
            log_msg(rs, "[NEWTON-RAPHSON] - Derivative too close to zero, aborting.");
            return ROOTSOLVER_NO_CONVERGENCE;
        }

        approx = approx_old - fx / dfx;

        err_rel = fabs(approx - approx_old) / (fabs(approx) + SAFETY_EPS) * 100.0;

        if (err_rel < tol)
        {
            log_msg(rs, "[NEWTON-RAPHSON] - Solution found: %.10f", approx);
            *sol = approx;
            return ROOTSOLVER_OK;
        }

        approx_old = approx;
        i += 1;
    }

    log_msg(rs, "[NEWTON-RAPHSON] - Method failed after %d iterations. The procedure was unsuccessful.", max_iter);
    return ROOTSOLVER_NO_CONVERGENCE;
}