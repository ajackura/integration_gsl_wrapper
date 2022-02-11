# integration_gsl_wrapper

**C++ wrapper for GSL numerical integration routines**

### requirements

- *c++20* for concepts type constraints
- *cmake 3.9* for compiling included tests
- *gnu scientific library 2.7* <https://www.gnu.org/software/gsl/>

### features

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