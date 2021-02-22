// #ifndef UNIFORM_HPP
// #define UNIFORM_HPP

// #include <Rcpp.h>
// using namespace Rcpp;

// class Uniform {
// public:
//   Uniform(double min_, double max_) :
//   min(min_), max(max_) {}
//   NumericVector draw(int n) const {
//     RNGScope scope;
//     return runif(n, min, max);
//   }
//   double min, max;
// };
// double uniformRange(Uniform* w) {
//   return w->max - w->min;
// }

// RCPP_MODULE(unif_module) {
//   class_<Uniform>("Uniform")
//   .constructor<double,double>()
//   .field("min", &Uniform::min)
//   .field("max", &Uniform::max)
//   .method("draw", &Uniform::draw)
//   .method("range", &uniformRange)
//   ;
// }

// #endif
