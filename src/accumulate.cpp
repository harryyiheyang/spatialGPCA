// SVGPC native core. Mathematical contract: inst/doc/MATHEMATICS.md.
#include "native.h"
#include <fstream>
void add_block(Context &c, mat Y) {
  auto t = Clock::now();
  Y.each_col() %= c.rootw;
  Y = residual(Y, c.Xq);
  c.energy += arma::accu(arma::square(Y));
  mat J = c.T.t() * Y;
  c.timing["projection"] += seconds(t);
  t = Clock::now();
  int n = J.n_rows, k = J.n_cols;
  double one = 1;
  char lower = 'L', nt = 'N';
  svgpc_dsyrk(&lower, &nt, &n, &k, &one, J.memptr(), &n, &one, c.cov.memptr(), &n);
  c.timing["covariance"] += seconds(t);
  c.p += Y.n_cols;
}
void finish(Context &c) {
  c.cov = arma::symmatl(c.cov);
  c.a = c.cov.diag();
  c.b = c.energy - arma::accu(c.a);
  if (c.b < -1e-9 * std::max(c.energy, 1.0))
    Rcpp::stop("Negative complement energy");
  c.b = c.rank() == c.r() ? 0.0 : std::max(0.0, c.b);
  c.ready = true;
}
void fresh(Context &c) {
  if (c.ready || c.p || c.read_variants)
    Rcpp::stop("Statistics already accumulated; create or load a fresh prepared model to avoid "
               "double counting");
}
// [[Rcpp::export]]
Rcpp::List cpp_accumulate_matrix(SEXP ptr, const arma::mat &G, int block_size) {
  auto c = context(ptr);
  fresh(*c);
  if (G.n_rows != (uword)c->n() || !G.n_cols || !G.is_finite() || block_size < 1)
    Rcpp::stop("Invalid location genotype matrix");
  for (uword j = 0; j < G.n_cols; j += block_size) {
    Rcpp::checkUserInterrupt();
    add_block(*c, G.cols(j, std::min(G.n_cols - 1, j + block_size - 1)));
  }
  finish(*c);
  return information(*c);
}

// [[Rcpp::export]]
Rcpp::List cpp_accumulate_file(SEXP ptr, std::string path, Rcpp::IntegerVector locations,
                               const arma::vec &af, int block_size) {
  auto c = context(ptr);
  fresh(*c);
  if (block_size < 1 || !af.n_elem || !af.is_finite() || arma::any(af <= 0) || arma::any(af >= 1))
    Rcpp::stop("Require positive block size and allele frequencies strictly between zero and one");
  vec counts(c->n(), arma::fill::zeros);
  for (int v : locations) {
    if (v < 0 || v > c->n())
      Rcpp::stop("Invalid sample to location mapping");
    if (v) counts[v - 1]++;
  }
  if (!arma::approx_equal(counts, c->weights, "absdiff", 0))
    Rcpp::stop("Location counts differ from the matched FAM samples");
  std::ifstream file(path, std::ios::binary);
  if (!file) Rcpp::stop("Cannot open BED file");
  unsigned char header[3];
  file.read(reinterpret_cast<char *>(header), 3);
  if (!file || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("Require a SNP-major PLINK BED file");
  const size_t ns = locations.size(), stride = (ns + 3) / 4;
  file.seekg(0, std::ios::end);
  if (file.tellg() != static_cast<std::streamoff>(3 + stride * af.n_elem))
    Rcpp::stop("BED size disagrees with FAM sample count or BIM variant count");
  file.seekg(3, std::ios::beg);
  std::vector<unsigned char> bytes(stride);
  vec scales = 1 / arma::sqrt(2 * af % (1 - af));
  for (uword start = 0; start < af.n_elem; start += block_size) {
    Rcpp::checkUserInterrupt();
    uword bs = std::min((uword)block_size, af.n_elem - start);
    mat G(c->n(), bs, arma::fill::zeros);
    auto t = Clock::now();
    for (uword j = 0; j < bs; ++j) {
      uword v = start + j;
      file.read(reinterpret_cast<char *>(bytes.data()), stride);
      if (!file) Rcpp::stop("Truncated BED at variant %d", static_cast<int>(v + 1));
      c->read_variants++;
      for (size_t i = 0; i < ns; ++i) {
        if (!locations[i]) continue;
        unsigned code = (bytes[i / 4] >> (2 * (i % 4))) & 3;
        double x;
        if (code == 1) {
          c->missing++;
          x = 2 * af[v];
        } else {
          x = code == 0 ? 2.0 : (code == 2 ? 1.0 : 0.0);
        }
        G(locations[i] - 1, j) += x;
      }
      G.col(j) = (G.col(j) / c->weights - 2 * af[v]) * scales[v];
    }
    c->timing["decode_aggregate_standardize"] += seconds(t);
    add_block(*c, std::move(G));
  }
  finish(*c);
  return information(*c);
}
