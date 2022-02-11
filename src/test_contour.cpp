#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>

#include "integration.hpp"

template <typename Tp>
std::complex<Tp> fcn(const std::complex<Tp> &z)
{
    return std::conj(z);
}

template <typename Tp>
std::complex<Tp> fcn2(const std::complex<Tp> &z)
{
    return Tp{1} / z;
}

template <typename Tp>
void test1(void)
{
    std::complex<Tp> i(Tp{0},Tp{1});
    std::vector<std::complex<Tp>> path{std::complex<Tp>(0),std::complex<Tp>(1),std::complex<Tp>(1)+i};
    auto ans = quad::integrate<Tp, decltype(fcn<Tp>)>(fcn<Tp>, path);
    std::cout << ans << std::endl;
}

template <typename Tp>
void test2(void)
{
    std::complex<Tp> i(Tp{0},Tp{1});
    std::vector<std::complex<Tp>> path{std::complex<Tp>(0),std::complex<Tp>(1)+i};
    auto ans = quad::integrate<Tp, decltype(fcn<Tp>)>(fcn<Tp>, path);
    std::cout << ans << std::endl;
}

template <typename Tp>
void test3(void)
{
    std::complex<Tp> i(Tp{0},Tp{1});
    std::vector<std::complex<Tp>> path{std::complex<Tp>(1),i,-std::complex<Tp>(1),-i,std::complex<Tp>(1)};
    auto ans = quad::integrate<Tp, decltype(fcn2<Tp>)>(fcn2<Tp>, path);
    std::cout << ans << std::endl;
}

template <typename Tp>
std::complex<Tp> func(Tp x)
{
    return std::exp(-std::complex<Tp>(x));
}

int main()
{
    test1<double>();
    test2<double>();
    test3<double>();

    return 0;
}