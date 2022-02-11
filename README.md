# integration_gsl_wrapper

**C++ wrapper for GSL numerical integration routines**

### requirements

- *c++20* for concepts type constraints
- *cmake 3.9* for compiling included tests
- *gnu scientific library 2.7* <https://www.gnu.org/software/gsl/>

### about

This is a simple single header file wrapper for the adaptive GSL integration routines. Provides and templated integrate function for any univariate function, which is overloaded for complex type, as well as a version to do contour integrals.

### header

*integration.hpp*
- contains a class which wraps the GSL initializations, as well as the primary user-callable functions to perform the integration.
`integrate(func, lower, upper, max_abs_err, max_rel_err, max_iter)`

- an overloaded verion is included to perform contour integrals, in which the user creates a vector of points *path* which specify the vertices of a straight-line path.
`integrate(func, path, max_abs_err, max_rel_err, max_iter)`

### tests

*test_integrals.hpp*
- few test cases

*test_contour.hpp*
- tests contour integrals


----

### compilation and building examples

Since this is a single header wrapper, simply include `integration.hpp` in your code and call the `integrate` functions. This wrapper was tested with g++11 on Mac OS X. The intel compiler as of the time of writing this, does not support *concepts*, which are necessary for type deduction for the numerical quadrature wrapper. To build the examples included, do the following:

1. clone repository, `cd` to root
2. run cmake: I do this from a script where I specify compilers

        mkdir build && cd build
        CC=/usr/local/bin/gcc-11 
        CXX=/usr/local/bin/g++-11 
        cmake ../

3. intsall tests: run `make && make install`, test executables are installed in `bin` in root.