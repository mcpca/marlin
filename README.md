marlin
===

A C++14 implementation of a Lax-Friedrichs fast sweeping method for first-order
static Hamilton-Jacobi equations.

The implementation is based on the algorithm described in Kao, Osher and Qian,
_Lax-Friedrichs sweeping scheme for static Hamilton-Jacobi equations_ (available
[here](https://doi.org/10.1016/j.jcp.2003.11.007)),
and is parallelized using the hyperplane stepping method described in Detrixhe,
Gibou and Min, _A parallel fast sweeping method for the Eikonal
equation_ ([link](https://doi.org/10.1016/j.jcp.2012.11.042)).

It is dimension-agnostic (the dimension is set at compile time).

Building
---

A C++14 compiler with OpenMP support is required to build the library.

```shell
  $ mkdir build && cd build
  $ cmake .. -DCMAKE_BUILD_TYPE=RELEASE
  $ make
```

Usage
---

The number of dimensions and the floating value type to use in computations
(i.e., `float`, `double` or `long double`) can be set in the header
`marlin/include/defs.hpp`.
The dimension can also be set using the macro `MARLIN_N_DIMS`.

The solver constructor takes four parameters:

- an array of floating point values which contains the values of the
  right-hand-side of the equation in row major order.
  Gridoints at which the Dirichlet boundary condition is imposed must have their
  value inverted in sign and shifted by negative one, e.g., if the boundary
  condition value is `0` at some point, the corresponding value in the array
  should be `-1`;

- an array containing the size of the grid in each dimension;

- the limits of the computational domain as an array of (min, max) pairs;

- the initial value of the solution and the tolerance (the epsilon value for the
  convergence criterion).
  These are passed inside a struct to avoid having a constructor which takes
  more than one argument of the same type, which can lead to errors in the order
  of arguments.

The Hamiltonian and the viscosity coefficients can be any kind of callable
object (function pointers, functors, lambdas, etc).

Check the examples for how to set up and call the solver.

Examples
---

Four examples are included in the `examples` directory.
You must have the [h5cpp](https://github.com/ess-dmsc/h5cpp) library installed
in order to build them.

In `examples/scripts` there are Python scripts which set up the problem, call
the solver and display the results.
For instance, to run the `eikonal3d` example:

```shell
  $ cd build
  $ make examples
  $ cd ../examples/scripts
  $ mkdir ../data
  $ ./eikonal3d.py
```

(`make examples` builds all the examples, so to run subsequent examples all you
have to do is call the corresponding script.)

Since the default number of dimensions is 3, you'll have to change the number of
dimensions to 2 and recompile in order to run the two-dimensional examples.

You'll need to have `numpy` and `h5py` installed to run the examples.
The two-dimensional examples also require `matplotlib`.

Acknowledgments
---

* [LSTS-FEUP](https://lsts.pt)
* Jorge Estrela da Silva
