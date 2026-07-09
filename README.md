# Celestial

Celestial is a numerical methods library written in C, built as a learning
project to study the language. **It's in a very early stage of
development**, so the API may change without notice.

## Current status

- **Root finding for nonlinear equations** (`roots_solver`)
  - Bisection
  - Regula falsi
  - Fixed point
  - Newton-Raphson
- **Initial value problems / ODEs** (`ivp_solver`)
  - Euler's method
  - 4th order Runge-Kutta (RK4)

Both modules support an optional logging callback (`LogFunc`) to report
errors or internal messages without depending on `printf`.

## Building

The project uses CMake and has no external dependencies beyond the C
standard library (and `libm` on non-Windows systems).

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

This produces a static library `celestial`. To use it in another CMake
project:

```cmake
add_subdirectory(celestial)
target_link_libraries(my_project PRIVATE celestial)
```

## Quick example

```c
#include <celestial/roots_solver.h>
#include <math.h>
#include <stdio.h>

double f(double x)  { return x * x - 2.0; }
double df(double x) { return 2.0 * x; }

int main(void)
{
    RootSolver rs = rootsolver_create(f, df, NULL);

    double sol;
    RootSolverStatus status = rootsolver_newtonraphson(&rs, &sol, 1.0, 1e-8, 100);

    if (status == ROOTSOLVER_OK)
        printf("Approximate root: %f\n", sol);

    return 0;
}
```

## Roadmap (ideas, no timeline commitment)

- Interpolation (Lagrange, splines)
- Numerical integration (trapezoidal rule, Simpson's rule)
- Linear systems of equations
- Curve fitting / least squares

## License

Not defined yet.
