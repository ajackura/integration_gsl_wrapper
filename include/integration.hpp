/*
 * Copyright (c) 2022 Andrew W. Jackura
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * Under Section 7 of GPL version 3, you are granted additional
 * permissions described in the GCC Runtime Library Exception, version
 * 3.1, as published by the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License and
 * a copy of the GCC Runtime Library Exception along with this program;
 * see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef QUAD_GSL_INTEGRATION_HPP__
#define QUAD_GSL_INTEGRATION_HPP__

/**
 * @file integration.hpp
 * @author Andrew W. Jackura (ajackura@jlab.org)
 * @brief wrapper for gsl qags integration routines implemented in <http://www.gnu.org/software/gsl/>.
 * currently, the adaptive routines supported are: QAGS, QAGIL, QAGIU, QAGI
 *
 * The base of this wrapper was taken from a stackexchange post
 * <https://scicomp.stackexchange.com/questions/20786/c-library-for-numerical-intergration-quadrature>
 * by user Henri Menke <https://scicomp.stackexchange.com/users/24680/henri-menke>.
 */

#include <cmath>
#include <complex>
#include <concepts>
#include <functional>
#include <memory>
#include <utility>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

namespace quad
{

    /**
     * @brief quadrature class GSL wrapper for adaptive integration
     */
    template <typename Tp, typename Func>
    class gsl_quad
    {
        Func func;
        std::size_t limit;
        std::unique_ptr<gsl_integration_workspace,
                        std::function<void(gsl_integration_workspace *)>>
            workspace;

        static Tp gsl_wrapper(Tp x, void *p)
        {
            gsl_quad *t = reinterpret_cast<gsl_quad *>(p);
            return t->func(x);
        }

    public:
        gsl_quad(Func func, std::size_t limit)
            : func(func), limit(limit), workspace(gsl_integration_workspace_alloc(limit), gsl_integration_workspace_free)
        {
        }

        Tp integrate(Tp min, Tp max, Tp epsabs, Tp epsrel)
        {
            gsl_set_error_handler_off();
            gsl_function gsl_f;
            gsl_f.function = &gsl_wrapper;
            gsl_f.params = this;

            Tp result, error;
            if (!std::isinf(min) && !std::isinf(max))
            {
                gsl_integration_qags(&gsl_f, min, max,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error);
            }
            else if (std::isinf(min) && !std::isinf(max))
            {
                gsl_integration_qagil(&gsl_f, max,
                                      epsabs, epsrel, limit,
                                      workspace.get(), &result, &error);
            }
            else if (!std::isinf(min) && std::isinf(max))
            {
                gsl_integration_qagiu(&gsl_f, min,
                                      epsabs, epsrel, limit,
                                      workspace.get(), &result, &error);
            }
            else
            {
                gsl_integration_qagi(&gsl_f,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error);
            }

            return result;
        }
    };

    // struct to test for complex types
    template <typename Tp>
    struct is_complex : std::false_type
    {
    };

    // struct to test for complex types
    template <std::floating_point Tp>
    struct is_complex<std::complex<Tp>> : std::true_type
    {
    };

    // complex concept
    template <typename Tp>
    concept complex_floating_point = is_complex<Tp>::value;

    /**
     * @brief integrator for floating point types
     *
     * @param func function of single variable to integrate
     * @param lower lower limit of integration
     * @param upper upper limit of integration
     * @param max_abs_err target absolute error
     * @param max_rel_err target relative error
     * @param max_iter max number of iterations
     * @return integral of func
     */
    template <typename Tp, typename Func>
    requires std::floating_point<std::invoke_result_t<Func, Tp>>
    auto integrate(const Func &func, 
                   const Tp &lower, 
                   const Tp &upper, 
                   const Tp &max_abs_err = Tp{1e-12}, 
                   const Tp &max_rel_err = Tp{1e-12}, 
                   const std::size_t &max_iter = 1024)
        -> std::invoke_result_t<Func, Tp>
    {
        return gsl_quad<Tp, Func>(func, max_iter).integrate(lower, upper, max_abs_err, max_rel_err);
    }

    /**
     * @brief integrator for complex floating point types, with real integration limits
     *
     * @param func function of single variable to integrate
     * @param lower lower limit of integration
     * @param upper upper limit of integration
     * @param max_abs_err target absolute error
     * @param max_rel_err target relative error
     * @param max_iter max number of iterations
     * @return integral of func
     */
    template <typename Tp, typename Func>
    requires complex_floating_point<std::invoke_result_t<Func, Tp>>
    std::complex<Tp>
    integrate(const Func &func,
              const Tp &lower, 
              const Tp &upper, 
              const Tp &max_abs_err = Tp{1e-12}, 
              const Tp &max_rel_err = Tp{1e-12}, 
              const std::size_t &max_iter = 1024)
    {
        const auto real = [&](Tp x) -> Tp
        { 
            return std::real(func(x)); 
        };
        const auto imag = [&](Tp x) -> Tp
        { 
            return std::imag(func(x)); 
        };

        const auto quad_real = integrate<Tp, decltype(real)>(real, lower, upper, max_abs_err, max_rel_err, max_iter);
        const auto quad_imag = integrate<Tp, decltype(imag)>(imag, lower, upper, max_abs_err, max_rel_err, max_iter);

        return {quad_real, quad_imag};
    }

    /**
     * @brief overloaded routine to do path integrals
     *
     * @param func function of single variable to integrate
     * @param path vector of points indicating vertices of path
     * @param max_abs_err target absolute error
     * @param max_rel_err target relative error
     * @param max_iter max number of iterations
     * @return integral of func
     */
    template <typename Tp, typename Func>
    requires complex_floating_point<std::invoke_result_t<Func, Tp>>
    std::complex<Tp>
    integrate(const Func &func,
              const std::vector<std::complex<Tp>> &path,
              const Tp &max_abs_err = Tp{1e-12},
              const Tp &max_rel_err = Tp{1e-12},
              const std::size_t &max_iter = 1024)
    {
        const auto lower = Tp{0};
        const auto upper = Tp{1};
        std::complex<Tp> result(Tp{0}, Tp{0});
        for (std::size_t i = 0; i < path.size() - 1; ++i)
        {
            const auto z0 = path[i];
            const auto dz = path[i + 1] - path[i];

            const auto func_trans = [&](Tp x) -> std::complex<Tp>
            {
                return func(z0 + x * dz) * dz;
            };
            result += integrate<Tp, decltype(func_trans)>(func_trans, lower, upper, max_abs_err, max_rel_err, max_iter);
        }
        return result;
    }

} // namespace quad

#endif // QUAD_GSL_INTEGRATION_HPP__