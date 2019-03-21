fsm
===

A C++14 implementation of a Lax-Friedrichs fast sweeping method for first-order
static Hamilton-Jacobi equations.

The implementation is based on the algorithm described in Kao, Osher and Qian,
_Lax-Friedrichs sweeping scheme for static Hamilton-Jacobi equations_ (available
[here](https://www.sciencedirect.com/science/article/pii/S0021999103006016)).

The implementation is dimension-agnostic (the dimension is set at compile time)
and uses the hdf5 format for data i/o.

Dependencies
---

* [h5cpp](https://github.com/ess-dmsc/h5cpp)

Building
---

```shell
  $ mkdir build && cd build
  $ cmake .. -DCMAKE_BUILD_TYPE=RELEASE
  $ make
```

Usage
---

The number of dimensions and the floating value type to use in computations
(i.e., `float`, `double` or `long double`) can be set in the header
`fsm/include/defs.hpp`.
The dimension can also be set using the macro `FSM_N_DIMS`.

The provided hdf5 file should have a dataset named `cost_function` with a number
of dimensions equal to that defined in `defs.hpp`.
This should contain an array of floating point values which give simultaneously
the values of the right-hand-side of the equation and the values at the points
on which the Dirichlet boundary condition is imposed.
This is accomplished by setting the boundary values to be negative and shifted
by negative one, e.g., if a boundary condition of `0` is to be set at some
point, the array in `cost_function` should have a value of `-1` at that point.

The Hamiltonian is passed to the solver as an object of any callable type (e.g.,
a function, function object, lambda or the result of a call to `std::bind`)
taking the index of a grid point and a `std::array` containing the values of the
derivative at that point.

Additionally, the grid 'vertices' (the limits of the computational boundary in
each dimension) and the artificial viscosity coefficients must be specified.

Check the examples for how to set up and call the solver.

Examples
---

Four examples are included in the `examples` directory.
These are Python scripts which set up the problem, call the solver and display
the results.
For instance, to run the `eikonal2d` example:

```shell
  $ cd build
  $ make examples
  $ cd ../examples/scripts
  $ python3 eikonal2d.py
```

(`make examples` builds all the examples, so to run subsequent examples all you
have to do is call the corresponding script.)

Since the default number of dimensions is 2, you'll have to change the number of
dimensions to 3 in order to run the `eikonal3d` example.

You'll need to have `numpy`, `h5py` and `matplotlib` installed to run the
examples.
