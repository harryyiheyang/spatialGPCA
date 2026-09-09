// SVGPC native core. Mathematical contract: inst/doc/MATHEMATICS.md.
#include "native.h"
vec shrink(const Context &c, double eta) {
  vec h(c.rank());
  for (int i = 0; i < c.rank(); ++i) {
    double x = eta - std::log(c.d[i]);
    h[i] = x > 0 ? std::exp(-x) / (1 + std::exp(-x)) : 1 / (1 + std::exp(x));
  }
  return h;
}
double criterion(const Context &c, double eta, const std::string &method) {
  if (!(c.energy > 0))
    Rcpp::stop("No residual variation");
  if (eta == INF)
    return 0;
  if (c.rank() == c.r()) {
    // Cancel powers of lambda analytically when no complement remains.
    // This handles the zero boundary without changing the independent-noise model.
    double lambda = std::exp(eta), energy = 0, inverse_trace = 0, logdet = 0;
    for (int i = 0; i < c.rank(); ++i) {
      double v = c.d[i] + lambda;
      energy += c.a[i] / v / (method == "GCV" ? v : 1.0);
      inverse_trace += 1 / v;
      logdet += std::log(v);
    }
    return method == "REML" ? logdet + c.r() * std::log(energy / c.energy)
                            : std::log(energy / c.energy) - 2 * std::log(inverse_trace / c.r());
  }
  double rss = c.b, logdet = 0, den = c.r() - c.rank();
  for (int i = 0; i < c.rank(); ++i) {
    double e = eta - std::log(c.d[i]);
    double z = e > 0 ? 1 / (1 + std::exp(-e)) : std::exp(e) / (1 + std::exp(e));
    den += z;
    rss += c.a[i] * (method == "REML" ? z : z * z);
    if (method == "REML") {
      double x = std::log(c.d[i]) - eta;
      logdet += std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x)));
    }
  }
  if (rss <= 0 || den <= 0)
    return INF;
  return method == "REML" ? logdet + c.r() * std::log(rss / c.energy)
                          : std::log(rss / c.energy) - 2 * std::log(den / c.r());
}
std::pair<double, double> bounded(const std::function<double(double)> &f, double a, double b) {
  const double g = (std::sqrt(5.0) - 1) / 2;
  double x = b - g * (b - a), y = a + g * (b - a), fx = f(x), fy = f(y);
  for (int it = 0; it < 150 && b - a > 1e-8; ++it) {
    if (fx < fy) {
      b = y;
      y = x;
      fy = fx;
      x = b - g * (b - a);
      fx = f(x);
    } else {
      a = x;
      x = y;
      fx = fy;
      y = a + g * (b - a);
      fy = f(y);
    }
  }
  return fx < fy ? std::make_pair(x, fx) : std::make_pair(y, fy);
}
// [[Rcpp::export]]
Rcpp::List cpp_fit(SEXP ptr, std::string method, Rcpp::CharacterVector operators, int components,
                   double fixed_lambda) {
  auto c = context(ptr);
  if (!c->ready)
    Rcpp::stop("Accumulate all variants before fitting");
  if (method != "REML" && method != "GCV")
    Rcpp::stop("Unknown selector");
  if (components < 1 || components > 50)
    Rcpp::stop("components must be 1..50");
  auto t = Clock::now();
  int calls = 0;
  auto objective = [&](double eta) {
    ++calls;
    return criterion(*c, eta, method);
  };
  double eta, best;
  std::string boundary;
  if (!std::isnan(fixed_lambda)) {
    if (fixed_lambda < 0)
      Rcpp::stop("lambda must be nonnegative");
    eta = std::log(fixed_lambda);
    best = objective(eta);
    boundary = "supplied";
  } else {
    double center = std::log(arma::median(c->d)), low = center - 30, step = .5;
    std::vector<double> val(121);
    for (int j = 0; j <= 120; ++j)
      val[j] = objective(low + j * step);
    eta = INF;
    best = 0;
    boundary = "infinity";
    if (method == "GCV" || c->rank() == c->r()) {
      double z = objective(-INF);
      if (z < best) {
        best = z;
        eta = -INF;
        boundary = "zero";
      }
    }
    for (int j = 1; j < 120; ++j)
      if (val[j] <= val[j - 1] && val[j] <= val[j + 1]) {
        auto opt = bounded(objective, low + (j - 1) * step, low + (j + 1) * step);
        if (opt.second < best - 1e-12) {
          best = opt.second;
          eta = opt.first;
          boundary = "finite";
        }
      }
    if (std::min(val.front(), val.back()) < best - 1e-8)
      Rcpp::stop("Lambda optimum reaches search boundary; extend bounds before accepting fit");
  }
  vec h = shrink(*c, eta);
  double opttime = seconds(t);
  Rcpp::List outputs;
  for (auto item : operators) {
    std::string op = Rcpp::as<std::string>(item);
    if (op != "raw_F_PCA" && op != "H2" && op != "H" && op != "2H-H2")
      Rcpp::stop("Unknown operator");
    if (!arma::any(h > 0)) {
      outputs[op] =
          op == "raw_F_PCA"
              ? Rcpp::List::create(Rcpp::_["pcs"] = mat(c->n(), 0), Rcpp::_["eigenvalues"] = vec())
              : Rcpp::List::create(Rcpp::_["weighted_scores"] = mat(c->n(), 0),
                                   Rcpp::_["raw_scores"] = mat(c->n(), 0));
      continue;
    }
    if (op == "raw_F_PCA") {
      mat RH = c->R.each_row() % h.t();
      mat covariance = RH * c->cov * RH.t() / (double)c->p;
      vec ev;
      mat U;
      eigen_top(covariance, components, ev, U);
      mat pcs = c->T * arma::solve(arma::trimatu(c->R), U);
      pcs.each_col() /= c->rootw;
      outputs[op] = Rcpp::List::create(Rcpp::_["pcs"] = pcs, Rcpp::_["eigenvalues"] = ev);
    } else {
      vec w = h;
      if (op == "H2")
        w = h % h;
      else if (op == "2H-H2")
        w = 2 * h - h % h;
      vec sw = arma::sqrt(w);
      mat covariance = c->cov;
      covariance.each_col() %= sw;
      covariance.each_row() %= sw.t();
      vec ev;
      mat U;
      eigen_top(covariance, components, ev, U);
      mat weighted = U.each_col() % sw;
      mat coord = c->cov * weighted;
      coord.each_col() %= h;
      coord.each_row() /= arma::sqrt(ev).t();
      mat scores = c->T * coord;
      mat raw = scores;
      raw.each_col() /= c->rootw;
      outputs[op] = Rcpp::List::create(Rcpp::_["weighted_scores"] = scores,
                                       Rcpp::_["raw_scores"] = raw);
    }
  }
  return Rcpp::List::create(
      Rcpp::_["lambda"] = std::exp(eta),
      Rcpp::_["selector"] = std::isnan(fixed_lambda) ? method : "fixed_lambda",
      Rcpp::_["boundary"] = boundary, Rcpp::_["objective_relative"] = best,
      Rcpp::_["evaluations"] = calls, Rcpp::_["optimization_seconds"] = opttime,
      Rcpp::_["total_seconds"] = seconds(t), Rcpp::_["spatial_df"] = arma::accu(h),
      Rcpp::_["outputs"] = outputs);
}
// [[Rcpp::export]]
arma::mat cpp_fitted(SEXP ptr, const arma::mat &G, double lambda, bool coefficients) {
  auto c = context(ptr);
  if (G.n_rows != (uword)c->n() || !G.is_finite() || std::isnan(lambda) || lambda < 0)
    Rcpp::stop("Invalid fitted block");
  mat J = c->T.t() * (G.each_col() % c->rootw);
  if (coefficients) {
    vec w = arma::sqrt(c->d) / (c->d + lambda);
    return c->N * c->V * (J.each_col() % w);
  }
  vec h = c->d / (c->d + lambda);
  mat fitted = c->T * (J.each_col() % h);
  fitted.each_col() /= c->rootw;
  return fitted;
}
