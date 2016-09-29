#include <boost/math/distributions/normal.hpp>
#include <random>

#include <omp.h>

using namespace arma;

double pnorm(double q, double a, double b) {
  boost::math::normal_distribution<double> distribution(a, b);
  return cdf(distribution, q);
}

vec pnorm_std(vec q) {
  int n = q.n_elem;
  vec p = zeros(n);
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    p(i) = pnorm(q(i), 0, 1);
  }
  return p;
}

vec log_lefttail_pnorm_std(vec q) {
  int n = q.n_elem;
  vec p = zeros(n);
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    p(i) = -q(i) * q(i) / 2 - log(abs(q(i)));
  }
  return p;
}

class C11RNG {
 public:
  std::mt19937 gen;
  std::vector<std::mt19937> gen_vec;
  C11RNG() {
    std::random_device std_rand;
    int tMax = omp_get_max_threads();
    gen_vec.resize(tMax);
    gen = std::mt19937(std_rand());
    for (int i = 0; i < tMax; ++i) {
      gen_vec[i] = std::mt19937(std_rand());
    }
  };

  double draw_gamma(double alpha, double beta, std::mt19937& gen) {
    std::gamma_distribution<> gamma1(alpha, 1. / beta);
    return gamma1(gen);
  }

  double draw_gamma(double alpha, double beta) {
    return draw_gamma(alpha, beta, gen);
  }

  double draw_chisq(double alpha, std::mt19937& gen) {
    std::gamma_distribution<> gamma1(alpha / 2, 2.0);
    return gamma1(gen);
  }

  double draw_chisq(double alpha) { return draw_chisq(alpha, gen); }

  vec draw_gamma(vec alpha, vec beta) {
    int N = alpha.n_elem;
    vec randvec(N);

#pragma omp parallel for
    for (int i = 0; i < N; i++) {
      int tid = omp_get_thread_num();
      randvec(i) = draw_gamma(alpha(i), beta(i), gen_vec[tid]);
    }
    return randvec;
  }

  int draw_binomial(int n, double p, std::mt19937& gen) {
    std::binomial_distribution<int> distribution(n, p);
    return distribution(gen);
  };

  int draw_binomial(int n, double p) { return draw_binomial(n, p, gen); }

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

  double rnorm(double m, double sd, std::mt19937& gen) {
    std::normal_distribution<double> distribution(m, sd);
    return distribution(gen);
  }

  double rnorm(double m, double sd) { return rnorm(m, sd, gen); }

  double rexp(double alpha, std::mt19937& gen) {
    std::exponential_distribution<> exp1(alpha);
    double r = exp1(gen);
    while (std::isnan(r)) {
      r = exp1(gen);
    }
    return r;
  }

  double rexp(double alpha) { return rexp(alpha, gen); }

  double runif(std::mt19937& gen) {
    std::uniform_real_distribution<> unif(0., 1.);
    return unif(gen);
  }

  double runif() { return runif(gen); }

  // truncated normal

  double rtruncnorm_std_lb(double lb, std::mt19937& gen) {
    double alpha = (lb + sqrt(lb * lb + 4)) / 2.0;

    double delta = 1;
    double u = 1.5;
    double z = 0;

    while (u > delta) {
      z = rexp(alpha, gen) + lb;

      if (lb < alpha) {
        delta = exp(-(alpha - z) * (alpha - z) / 2);
      } else {
        delta = exp(-(alpha - z) * (alpha - z) / 2 +
                    (lb - alpha) * (lb - alpha) / 2);
      }
      u = runif(gen);
    }
    return z;
  }

  double rtruncnorm_std_ub(double ub, std::mt19937& gen) {
    return -rtruncnorm_std_lb(-ub, gen);
  }

  double rtruncnorm_std_lb_ub(double lb, double ub, std::mt19937& gen) {
    double trunc = 100;
    if (ub > trunc) ub = trunc;
    if (lb < -trunc) lb = -trunc;

    double z = 0;

    if (lb > ub) {
      z = lb;
    } else if (lb < trunc && ub > -trunc) {
      double delta = 1;
      double u = 1.5;

      while (u > delta) {
        z = runif(gen) * (ub - lb) + lb;

        if (ub < 0) {
          delta = exp((ub * ub - z * z) / 2);
        }
        if (lb > 0) {
          delta = exp((lb * lb - z * z) / 2);
        }
        if (lb <= 0 && ub >= 0) {
          delta = exp(-z * z / 2);
        }

        u = runif(gen);
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

  double rtruncnorm_std(double lb, double ub, std::mt19937& gen) {
    double r = 0;
    if (lb == -INFINITY && ub == INFINITY) {
      r = rnorm(0, 1, gen);
    } else if (lb == -INFINITY && ub < INFINITY) {
      r = rtruncnorm_std_ub(ub, gen);
    } else if (ub == INFINITY && lb > -INFINITY) {
      r = rtruncnorm_std_lb(lb, gen);
    } else {
      r = rtruncnorm_std_lb_ub(lb, ub, gen);
    }
    return r;
  }

  double rtruncnorm(double mean, double sigma, double lb, double ub,
                    std::mt19937& gen) {
    double lb1 = (lb - mean) / sigma;
    double ub1 = (ub - mean) / sigma;

    double r = rtruncnorm_std(lb1, ub1, gen);
    r = r * sigma + mean;

    while (std::isnan(r) || r < lb || r > ub) {
      r = rtruncnorm(mean, sigma, lb, ub, gen);
    }

    return r;
  }

  vec rtruncnorm(vec mean, vec sigma, vec lb, vec ub) {
    int n = mean.n_elem;
    vec r = zeros(n);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      int tid = omp_get_thread_num();
      r(i) = rtruncnorm(mean(i), sigma(i), lb(i), ub(i), gen_vec[tid]);
    }
    return r;
  }
};