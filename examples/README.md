Examples
===

This is a collection of four programs meant to illustrate the usage of the
solver.

eikonal2d
---

The simplest example of a static Hamilton-Jacobi equation.
A zero boundary condition is set at some grid point.
The equation is solved over the square `[-1, 1]x[-1, 1]`, with an adjustable
resolution.

Since the analytical solution to this equation is known, the maximum error among
all gridpoints is printed after the solver is called.

eikonal3d
---

The simplest example of a static Hamilton-Jacobi equation.
A zero boundary condition is set at some grid point.
The equation is solved over the square `[-1, 1]x[-1, 1]x[-1, 1]`, with an
adjustable resolution.

Since the analytical solution to this equation is known, the maximum error among
all gridpoints is printed after the solver is called.

Remember to compile the library with `marlin::dim` set to `3` in `defs.hpp` before
running this example.

**NOTE**: This example is meant to illustrate the library's capability to solve
problems in an arbitrary number of dimensions.
Since any three-dimensional problem involves a relatively high number of
gridpoints, this example may take a while to run.

constant_field
---

Computes the value function associated to a minimum time control problem for a
holonomic vehicle moving in two dimensions with a constant velocity disturbance.

The components of the disturbance and the target point are adjustable.
For zero disturbance the problem reduces to an eikonal equation.

geodesic
---

Computes the minimum time from a point to every other point on a surface.
This example was taken from Sethian and Vladimirsky, _Ordered Upwind Methods for
Static Hamilton-Jacobi Equations_ (available
[here](https://epubs.siam.org/doi/abs/10.1137/S0036142901392742)).
