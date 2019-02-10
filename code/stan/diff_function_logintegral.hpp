#include <boost/math/quadrature/exp_sinh.hpp>

template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
  diff_function_logintegral(const T0__& alpha,
                          const T1__& z,
                          const T2__& w,
                          const T3__& m, std::ostream* pstream__){
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type T;
    double inf = std::numeric_limits<double>::infinity();
    boost::math::quadrature::exp_sinh<double> integrator;
    auto f = [&](double x) {
      return std::pow(x, -value_of(alpha))  *  std::pow(value_of(z), x-1) * std::pow(value_of(w), std::pow(x-1, 2)); 
    };
    double termination = std::pow(std::numeric_limits<double>::epsilon(), 0.5);
    double error = 0.0;
    double L1 = 0.0;
    size_t levels;
    double Q = integrator.integrate(f, value_of(m), inf, termination, &error, &L1, &levels);
    T ans = std::log(Q); 
    return ans;
}