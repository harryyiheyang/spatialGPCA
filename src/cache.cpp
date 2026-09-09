// SVGPC native core. Mathematical contract: inst/doc/MATHEMATICS.md.
#include "native.h"
Rcpp::XPtr<Context> context(SEXP x) {
  Rcpp::XPtr<Context> p(x);
  if (!p)
    Rcpp::stop("Invalid context; use svgpc_load() to reload saved models");
  return p;
}
Rcpp::List information(Context &c) {
  return Rcpp::List::create(Rcpp::_["n_locations"] = c.n(), Rcpp::_["knots"] = c.knots.n_rows,
                            Rcpp::_["fixed_rank"] = c.Xq.n_cols, Rcpp::_["spatial_rank"] = c.rank(),
                            Rcpp::_["rho"] = c.rho,
                            Rcpp::_["n_snps"] = (double)c.p,
                            Rcpp::_["decoded_variants"] = (double)c.read_variants,
                            Rcpp::_["missing_genotypes"] = (double)c.missing,
                            Rcpp::_["statistics_ready"] = c.ready, Rcpp::_["seconds"] = c.timing);
}
// [[Rcpp::export]]
Rcpp::List cpp_info(SEXP ptr) {
  auto c = context(ptr);
  return information(*c);
}
// [[Rcpp::export]]
Rcpp::List cpp_statistics(SEXP ptr) {
  auto c = context(ptr);
  return Rcpp::List::create(Rcpp::_["covariance_sum"] = c->cov, Rcpp::_["total_energy"] = c->energy,
                            Rcpp::_["b"] = c->b, Rcpp::_["p"] = (double)c->p, Rcpp::_["d"] = c->d,
                            Rcpp::_["residual_df"] = c->r());
}
// [[Rcpp::export]]
Rcpp::List cpp_snapshot(SEXP ptr) {
  auto c = context(ptr);
  Rcpp::List out;
  out["version"] = 3;
#define SAVE(name) out[#name] = c->name;
  SAVE(xy)
  SAVE(knots) SAVE(C) SAVE(Xq) SAVE(N) SAVE(V) SAVE(T) SAVE(R) SAVE(cov) SAVE(weights)
      SAVE(d) SAVE(rho) SAVE(energy) SAVE(ready) SAVE(timing)
#undef SAVE
          out["p"] = (double)c->p;
  out["missing"] = (double)c->missing;
  out["read_variants"] = (double)c->read_variants;
  return out;
}
// [[Rcpp::export]]
SEXP cpp_restore(Rcpp::List s) {
  int version = Rcpp::as<int>(s["version"]);
  if (version != 2 && version != 3)
    Rcpp::stop("Unknown cache version");
  std::unique_ptr<Context> c(new Context);
#define LOADMAT(name) c->name = Rcpp::as<mat>(s[#name]);
  LOADMAT(xy)
  LOADMAT(knots) LOADMAT(C) LOADMAT(Xq) LOADMAT(N) LOADMAT(V) LOADMAT(T) LOADMAT(R)
      LOADMAT(cov)
#undef LOADMAT
          c->weights = Rcpp::as<vec>(s["weights"]);
  c->d = Rcpp::as<vec>(s["d"]);
  c->rootw = arma::sqrt(c->weights);
  c->rho = Rcpp::as<double>(s["rho"]);
  c->energy = Rcpp::as<double>(s["energy"]);
  c->p = Rcpp::as<double>(s["p"]);
  c->ready = Rcpp::as<bool>(s["ready"]);
  Rcpp::NumericVector tm = s["timing"];
  Rcpp::CharacterVector tn = tm.names();
  for (int i = 0; i < tm.size(); ++i)
    c->timing[Rcpp::as<std::string>(tn[i])] = tm[i];
  c->missing = Rcpp::as<double>(s["missing"]);
  c->read_variants = Rcpp::as<double>(s["read_variants"]);
  if (c->ready)
    finish(*c);
  return Rcpp::XPtr<Context>(c.release(), true);
}
