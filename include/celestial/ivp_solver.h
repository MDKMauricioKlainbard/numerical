/**
 * @file ivp_solver.h
 * @brief Numerical integrators for initial value problems (systems of ODEs).
 *
 * Solves systems of the form dy/dt = f(t, y) using single-step methods.
 * An @ref IVPSolver holds the scratch buffers needed by the higher-order
 * methods so that stepping functions can run without allocating memory.
 */

#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H

/**
 * @brief Signature for the right-hand side of the ODE system dy/dt = f(t, y).
 * @param t    Current independent variable (time).
 * @param y    Current state vector, of length @p n.
 * @param dydt Output vector, of length @p n, to be filled with dy/dt at (t, y).
 * @param n    Dimension of the system.
 */
typedef void (*OdeFunc)(double t, const double *y, double *dydt, int n);

/**
 * @brief Signature for the optional logging callback.
 * @param msg Null-terminated message describing progress or an error.
 */
typedef void (*LogFunc)(const char *msg);

/**
 * @brief Context and scratch space for stepping an ODE system.
 *
 * Create with @ref ivpsolver_create and release with @ref ivpsolver_destroy.
 * The `k1`..`k4` and `y_temp` buffers are internal working memory sized for
 * an `n`-dimensional system and are reused across every call to
 * @ref ivpsolver_step_euler / @ref ivpsolver_step_rk4.
 */
typedef struct
{
    OdeFunc f;      /**< Right-hand side of the ODE system, dy/dt = f(t, y). */
    int n;           /**< Dimension of the system (length of the state vector). */
    LogFunc log;     /**< Optional logging callback. May be NULL to disable logging. */
    double *k1;      /**< Internal scratch buffer (RK4 stage 1). */
    double *k2;      /**< Internal scratch buffer (RK4 stage 2). */
    double *k3;      /**< Internal scratch buffer (RK4 stage 3). */
    double *k4;      /**< Internal scratch buffer (RK4 stage 4). */
    double *y_temp;  /**< Internal scratch buffer for intermediate state evaluations. */
} IVPSolver;

/**
 * @brief Result/status code returned by @ref ivpsolver_create.
 */
typedef enum
{
    IVPSOLVER_OK = 0,           /**< The solver was initialized successfully. */
    IVPSOLVER_ALLOC_FAILED = 1  /**< Allocation of an internal scratch buffer failed. */
} IVPSolverStatus;

/**
 * @brief Initialize an IVP solver, allocating its internal scratch buffers.
 * @param ivp Solver instance to initialize.
 * @param f   Right-hand side of the ODE system. Must not be NULL.
 * @param n   Dimension of the system (number of state variables).
 * @param log Optional logging callback, or NULL to disable logging.
 * @return ::IVPSOLVER_OK on success, or ::IVPSOLVER_ALLOC_FAILED if a scratch
 *         buffer could not be allocated.
 * @note On failure, call @ref ivpsolver_destroy to safely release any
 *       buffers that were allocated before the failure.
 */
IVPSolverStatus ivpsolver_create(IVPSolver *ivp, OdeFunc f, int n, LogFunc log);

/**
 * @brief Advance the state vector by one step using the explicit Euler method.
 *
 * First-order accurate; intended for quick experiments rather than precise
 * integration. @p y is updated in place.
 *
 * @param ivp Initialized solver instance.
 * @param t   Current time.
 * @param y   State vector of length `ivp->n`, updated in place to y(t + dt).
 * @param dt  Step size.
 *
 * @par Example
 * Integrate dy/dt = -y (exponential decay) from t = 0 to t = 1:
 * @code
 * void decay(double t, const double *y, double *dydt, int n)
 * {
 *     (void)t; (void)n;
 *     dydt[0] = -y[0];
 * }
 *
 * IVPSolver ivp;
 * ivpsolver_create(&ivp, decay, 1, NULL);
 *
 * double y[1] = { 1.0 };
 * double dt = 0.01;
 *
 * for (double t = 0.0; t < 1.0; t += dt)
 *     ivpsolver_step_euler(&ivp, t, y, dt);
 *
 * ivpsolver_destroy(&ivp);
 * @endcode
 *
 * @warning Euler's method is only first-order accurate and can become
 *          numerically unstable for stiff systems or step sizes @p dt that
 *          are too large relative to the system's dynamics. Prefer
 *          @ref ivpsolver_step_rk4 unless raw speed matters more than
 *          accuracy.
 * @warning @p y must point to a buffer of exactly `ivp->n` doubles, matching
 *          the dimension passed to @ref ivpsolver_create. A shorter buffer
 *          results in a buffer overrun.
 */
void ivpsolver_step_euler(const IVPSolver *ivp, double t, double *y, double dt);

/**
 * @brief Advance the state vector by one step using classic 4th-order Runge-Kutta (RK4).
 *
 * Fourth-order accurate. @p y is updated in place.
 *
 * @param ivp Initialized solver instance.
 * @param t   Current time.
 * @param y   State vector of length `ivp->n`, updated in place to y(t + dt).
 * @param dt  Step size.
 *
 * @par Example
 * Integrate a 2D harmonic oscillator (y[0] = position, y[1] = velocity)
 * from t = 0 to t = 10:
 * @code
 * void oscillator(double t, const double *y, double *dydt, int n)
 * {
 *     (void)t; (void)n;
 *     dydt[0] = y[1];
 *     dydt[1] = -y[0];
 * }
 *
 * IVPSolver ivp;
 * if (ivpsolver_create(&ivp, oscillator, 2, NULL) != IVPSOLVER_OK)
 *     return;
 *
 * double y[2] = { 1.0, 0.0 }; // position = 1, velocity = 0
 * double dt = 0.01;
 *
 * for (double t = 0.0; t < 10.0; t += dt)
 *     ivpsolver_step_rk4(&ivp, t, y, dt);
 *
 * ivpsolver_destroy(&ivp);
 * @endcode
 *
 * @warning @p y must point to a buffer of exactly `ivp->n` doubles, matching
 *          the dimension passed to @ref ivpsolver_create. A shorter buffer
 *          results in a buffer overrun.
 * @warning Not reentrant/thread-safe for a shared @ref IVPSolver instance:
 *          each step reuses the solver's internal `k1`..`k4` and `y_temp`
 *          scratch buffers, so the same @p ivp must not be stepped
 *          concurrently from multiple threads. Use a separate @ref IVPSolver
 *          per thread if you need to integrate systems in parallel.
 */
void ivpsolver_step_rk4(const IVPSolver *ivp, double t, double *y, double dt);

/**
 * @brief Release the internal scratch buffers held by an IVP solver.
 * @param ivp Solver instance to tear down. Safe to call even if some
 *            buffers failed to allocate in @ref ivpsolver_create.
 */
void ivpsolver_destroy(IVPSolver *ivp);

#endif
