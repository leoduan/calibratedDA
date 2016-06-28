// #include <boost/math/distributions/normal.hpp>
#include <random>

using namespace arma;

// double pnorm(double q, double a, double b) {
//   boost::math::normal_distribution<double> distribution(a, b);
//   return cdf(distribution, q);
// }

// vec pnorm_std(vec q) {
//   int n = q.n_elem;
//   vec p = zeros(n);
//   for (int i = 0; i < n; ++i) {
//     p(i) = pnorm(q(i), 0, 1);
//   }
//   return p;
// }

class C11RNG {
 public:
  std::mt19937 gen;
  C11RNG() {
    std::random_device std_rand;
    gen = std::mt19937(std_rand());
  };

  int draw_binomial(int n, double p) {
    std::binomial_distribution<int> distribution(n, p);
    return distribution(gen);
  };

  arma::vec draw_binomial(int n, double p, int size) {
    std::binomial_distribution<int> distribution(n, p);

    arma::vec result(size);
    for (int i = 0; i < size; ++i) {
      result(i) = (double)distribution(gen);
    }
    return result;
  };

  arma::vec draw_binomial(int n, arma::vec p, int size) {
    arma::vec result(size);
    for (int i = 0; i < size; ++i) {
      if (p(i) < 1) {
        std::binomial_distribution<int> distribution(n, p(i));
        result(i) = (double)distribution(gen);
      } else {
        result(i) = n;
      }
    }
    return result;
  };

  double rnorm(double m, double sd) {
    std::normal_distribution<double> distribution(m, sd);
    return distribution(gen);
  }

  double rexp(double alpha) {
    std::exponential_distribution<> exp1(alpha);
    double r = exp1(gen);
    while (std::isnan(r)) {
      r = exp1(gen);
    }
    return r;
  }

  double runif() {
    std::uniform_real_distribution<> unif(0., 1.);
    return unif(gen);
  }

  // truncated normal

  double rtruncnorm_std_lb(double lb) {
    double alpha = (lb + sqrt(lb * lb + 4)) / 2.0;

    double delta = 1;
    double u = 1.5;
    double z = 0;

    while (u > delta) {
      z = rexp(alpha) + lb;

      if (lb < alpha) {
        delta = exp(-(alpha - z) * (alpha - z) / 2);
      } else {
        delta = exp(-(alpha - z) * (alpha - z) / 2 +
                    (lb - alpha) * (lb - alpha) / 2);
      }
      u = runif();
    }
    return z;
  }

  double rtruncnorm_std_ub(double ub) { return -rtruncnorm_std_lb(-ub); }

  double rtruncnorm_std_lb_ub(double lb, double ub) {
    double trunc = 6;
    if (ub > trunc) ub = trunc;
    if (lb < -trunc) lb = -trunc;

    double z = 0;

    if (lb > ub) {
      z = lb;
    } else if (lb < trunc && ub > -trunc) {
      double delta = 1;
      double u = 1.5;

      while (u > delta) {
        z = runif() * (ub - lb) + lb;

        if (ub < 0) {
          delta = exp((ub * ub - z * z) / 2);
        }
        if (lb > 0) {
          delta = exp((lb * lb - z * z) / 2);
        }
        if (lb <= 0 && ub >= 0) {
          delta = exp(-z * z / 2);
        }

        u = runif();
      }
    } else {
      if (fabs(lb) < fabs(ub)) {
        z = lb;
      } else {
        z = ub;
      }
    }
    return z;
  }

  double rtruncnorm_std(double lb, double ub) {
    double r = 0;
    if (lb == -INFINITY && ub == INFINITY) {
      r = rnorm(0, 1);
    } else if (lb == -INFINITY && ub < INFINITY) {
      r = rtruncnorm_std_ub(ub);
    } else if (ub == INFINITY && lb > -INFINITY) {
      r = rtruncnorm_std_lb(lb);
    } else {
      r = rtruncnorm_std_lb_ub(lb, ub);
    }
    return r;
  }

  double rtruncnorm(double mean, double sigma, double lb, double ub) {
    double lb1 = (lb - mean) / sigma;
    double ub1 = (ub - mean) / sigma;

    double r = rtruncnorm_std(lb1, ub1);
    r = r * sigma + mean;

    while (std::isnan(r) || r < lb || r > ub) {
      r = rtruncnorm(mean, sigma, lb, ub);
    }

    return r;
  }

  vec rtruncnorm(vec mean, vec sigma, vec lb, vec ub) {
    int n = mean.n_elem;
    vec r = zeros(n);
    for (int i = 0; i < n; ++i) {
      r(i) = rtruncnorm(mean(i), sigma(i), lb(i), ub(i));
    }
    return r;
  }
};