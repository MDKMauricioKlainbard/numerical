/**
 * @file roots_solver.h
 * @brief Root-finding methods for scalar nonlinear equations f(x) = 0.
 *
 * Provides bracketing (closed) methods -- bisection and regula falsi -- and
 * iterative (open) methods -- fixed point and Newton-Raphson. All methods
 * share the same @ref RootSolver context and report progress/errors through
 * an optional logging callback instead of writing to stdout/stderr directly.
 */

#ifndef ROOTS_SOLVER_H
#define ROOTS_SOLVER_H

/**
 * @brief Signature for the function whose root is being sought.
 *
 * Evaluated as `f(x)`, where `x` is the point at which to evaluate the function.
 *
 * @return f(x).
 */
typedef double (*FuncPtr)(double);

/**
 * @brief Signature for the optional logging callback.
 * @param msg Null-terminated message describing progress or an error.
 */
typedef void (*LogFunc)(const char *msg);

/**
 * @brief Context shared by every root-finding method.
 *
 * Create one with @ref rootsolver_create and pass its address to any of the
 * `rootsolver_*` functions.
 */
typedef struct
{
    FuncPtr f;   /**< Function f(x) whose root is being sought. Required. */
    FuncPtr df;  /**< Analytical derivative f'(x). May be NULL; when NULL,
                      methods that need it (Newton-Raphson) fall back to a
                      finite-difference approximation. */
    LogFunc log; /**< Optional logging callback. May be NULL to disable logging. */
} RootSolver;

/**
 * @brief Result/status code returned by every `rootsolver_*` solving function.
 */
typedef enum
{
    ROOTSOLVER_OK = 0,              /**< A root was found within tolerance. */
    ROOTSOLVER_INVALID_INTERVAL = 1,/**< The provided bracket [a, b] is invalid or does not bracket a root. */
    ROOTSOLVER_NO_CONVERGENCE = 2   /**< The method did not converge within `max_iter` iterations. */
} RootSolverStatus;

/**
 * @brief Build a new root-solver context.
 * @param f   Function whose root is sought. Must not be NULL.
 * @param df  Analytical derivative of @p f, or NULL to use a numerical
 *            approximation where needed.
 * @param log Optional logging callback, or NULL to disable logging.
 * @return An initialized @ref RootSolver, by value.
 */
RootSolver rootsolver_create(FuncPtr f, FuncPtr df, LogFunc log);

/* Closed (bracketing) methods */

/**
 * @brief Find a root of f in [a, b] using the bisection method.
 *
 * Requires f(a) and f(b) to have opposite signs.
 *
 * @param rs       Root-solver context.
 * @param sol      Output parameter receiving the approximated root on success.
 * @param a        Lower bound of the search interval.
 * @param b        Upper bound of the search interval.
 * @param tol      Relative error tolerance (percentage) used as the stopping criterion.
 * @param max_iter Maximum number of iterations allowed.
 * @return ::ROOTSOLVER_OK on success, ::ROOTSOLVER_INVALID_INTERVAL if the
 *         interval does not bracket a root, or ::ROOTSOLVER_NO_CONVERGENCE
 *         if the iteration limit is reached.
 */
RootSolverStatus rootsolver_bisection(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter);

/**
 * @brief Find a root of f in [a, b] using the regula falsi (false position) method.
 *
 * Requires f(a) and f(b) to have opposite signs. Uses the Illinois-style
 * halving of a stale endpoint's function value to avoid stagnation.
 *
 * @param rs       Root-solver context.
 * @param sol      Output parameter receiving the approximated root on success.
 * @param a        Lower bound of the search interval.
 * @param b        Upper bound of the search interval.
 * @param tol      Relative error tolerance (percentage) used as the stopping criterion.
 * @param max_iter Maximum number of iterations allowed.
 * @return ::ROOTSOLVER_OK on success, ::ROOTSOLVER_INVALID_INTERVAL if the
 *         interval does not bracket a root, or ::ROOTSOLVER_NO_CONVERGENCE
 *         if the iteration limit is reached.
 */
RootSolverStatus rootsolver_regulafalsi(const RootSolver *rs, double *const sol, double a, double b, const double tol, const int max_iter);

/* Open (iterative) methods */

/**
 * @brief Find a fixed point of f (i.e. x such that f(x) = x) by direct iteration.
 *
 * @param rs       Root-solver context.
 * @param sol      Output parameter receiving the approximated fixed point on success.
 * @param approx   Initial guess.
 * @param tol      Relative error tolerance (percentage) used as the stopping criterion.
 * @param max_iter Maximum number of iterations allowed.
 * @return ::ROOTSOLVER_OK on success, or ::ROOTSOLVER_NO_CONVERGENCE if the
 *         iteration limit is reached.
 */
RootSolverStatus rootsolver_fixedpoint(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter);

/**
 * @brief Find a root of f starting from an initial guess using Newton-Raphson.
 *
 * Uses @ref RootSolver::df when available; otherwise approximates the
 * derivative numerically via a forward finite difference.
 *
 * @param rs       Root-solver context.
 * @param sol      Output parameter receiving the approximated root on success.
 * @param approx   Initial guess.
 * @param tol      Relative error tolerance (percentage) used as the stopping criterion.
 * @param max_iter Maximum number of iterations allowed.
 * @return ::ROOTSOLVER_OK on success, or ::ROOTSOLVER_NO_CONVERGENCE if the
 *         derivative becomes numerically zero or the iteration limit is reached.
 */
RootSolverStatus rootsolver_newtonraphson(const RootSolver *rs, double *const sol, double approx, const double tol, const int max_iter);
#endif
