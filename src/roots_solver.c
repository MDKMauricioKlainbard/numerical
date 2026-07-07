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

static RootSolverStatus validate_bracket(const RootSolver *rs, double a, double b, double *fa, double *fb, const char *method_name)
{
    if (a >= b)
    {
        log_msg(rs, "[%s] - Invalid interval: a must be less than b", method_name);
        return ROOTSOLVER_INVALID_INTERVAL;
    }

    *fa = rs->f(a);
    *fb = rs->f(b);

    if ((*fa * *fb) > 0.0)
    {
        log_msg(rs, "[%s] - No solution in interval [%.10f, %.10f]", method_name, a, b);
        return ROOTSOLVER_INVALID_INTERVAL;
    }

    return ROOTSOLVER_OK;
}

static RootSolverStatus report_success(const RootSolver *rs, double *const sol, double approx, const char *method_name)
{
    log_msg(rs, "[%s] - Solution found: %.10f", method_name, approx);
    *sol = approx;
    return ROOTSOLVER_OK;
}

static RootSolverStatus report_failure(const RootSolver *rs, const char *method_name, int max_iter)
{
    log_msg(rs, "[%s] - Method failed after %d iterations. The procedure was unsuccessful.", method_name, max_iter);
    return ROOTSOLVER_NO_CONVERGENCE;
}

static double relative_error(double approx, double approx_old)
{
    return fabs(approx - approx_old) / (fabs(approx) + SAFETY_EPS) * 100.0;
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
    double fa, fb;
    RootSolverStatus status = validate_bracket(rs, a, b, &fa, &fb, "BISECTION");
    if (status != ROOTSOLVER_OK)
        return status;

    int i = 1;
    double approx_old = a;

    double approx, fapprox;

    while (i <= max_iter)
    {
        approx = (a + b) / 2.0;
        fapprox = rs->f(approx);

        if (relative_error(approx, approx_old) < tol)
            return report_success(rs, sol, approx, "BISECTION");

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

    return report_failure(rs, "BISECTION", max_iter);
}

RootSolverStatus rootsolver_regulafalsi(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter)
{
    double fa, fb;

    RootSolverStatus status = validate_bracket(rs, a, b, &fa, &fb, "REGULA FALSI");
    if (status != ROOTSOLVER_OK)
        return status;

    int i = 1;
    double approx_old = a;

    double approx, fapprox;

    int left_repeat = 0;
    int right_repeat = 0;

    while (i <= max_iter)
    {
        approx = b - (fb * (a - b)) / (fa - fb);
        fapprox = rs->f(approx);

        if (relative_error(approx, approx_old) < tol)
            return report_success(rs, sol, approx, "REGULA FALSI");

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

    return report_failure(rs, "REGULA FALSI", max_iter);
}

RootSolverStatus rootsolver_fixedpoint(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter)
{
    double approx_old = approx;
    int i = 1;

    while (i <= max_iter)
    {
        approx = rs->f(approx_old);

        if (relative_error(approx, approx_old) < tol)
            return report_success(rs, sol, approx, "FIXED POINT");

        approx_old = approx;
        i += 1;
    }

    return report_failure(rs, "FIXED POINT", max_iter);
}

RootSolverStatus rootsolver_newtonraphson(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter)
{
    int i = 1;
    double approx_old = approx;

    double fx, dfx;

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

        if (relative_error(approx, approx_old) < tol)
            return report_success(rs, sol, approx, "NEWTON-RAPHSON");

        approx_old = approx;
        i += 1;
    }

    return report_failure(rs, "NEWTON-RAPHSON", max_iter);
}