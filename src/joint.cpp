// Individual-level BED block decoder for the joint U + SPDE decomposition.
// Returns standardized genotypes (g - 2f) / sqrt(2f(1-f)); missing calls become 0 (= 2f).
#include <RcppArmadillo.h>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat cpp_joint_bed_block(std::string path, int file_samples, Rcpp::IntegerVector rows,
                              int start, int count, const arma::vec &af, int threads) {
  if (file_samples < 1 || start < 0 || count < 1 || (arma::uword)(start + count) > af.n_elem)
    Rcpp::stop("Invalid BED block request");
  const size_t stride = ((size_t)file_samples + 3) / 4;
  std::ifstream file(path, std::ios::binary);
  if (!file) Rcpp::stop("Cannot open BED file");
  unsigned char header[3];
  file.read(reinterpret_cast<char *>(header), 3);
  if (!file || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("Require a SNP-major PLINK BED file");
  std::vector<unsigned char> bytes(stride * (size_t)count);
  file.seekg(3 + (std::streamoff)stride * start, std::ios::beg);
  file.read(reinterpret_cast<char *>(bytes.data()), (std::streamsize)bytes.size());
  if (!file) Rcpp::stop("Truncated BED");
  const int N = rows.size();
  std::vector<int> r(N);
  for (int i = 0; i < N; ++i) {
    if (rows[i] < 1 || rows[i] > file_samples) Rcpp::stop("Invalid FAM row index");
    r[i] = rows[i] - 1;
  }
  arma::mat X(N, count);
#ifdef _OPENMP
  if (threads < 1) threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < count; ++j) {
    const unsigned char *b = bytes.data() + stride * (size_t)j;
    const double f = af[start + j], mu = 2 * f, sc = 1 / std::sqrt(2 * f * (1 - f));
    // PLINK codes: 00 hom A1 (2 copies), 01 missing, 10 het, 11 hom A2.
    const double lut[4] = {(2 - mu) * sc, 0.0, (1 - mu) * sc, (0 - mu) * sc};
    double *col = X.colptr(j);
    for (int i = 0; i < N; ++i) {
      const int k = r[i];
      col[i] = lut[(b[k >> 2] >> (2 * (k & 3))) & 3];
    }
  }
  return X;
}

// [[Rcpp::export]]
int cpp_joint_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// Location sums Z'X and total sum of squares of an N x b block in one pass.
// [[Rcpp::export]]
Rcpp::List cpp_joint_aggregate(const arma::mat &X, Rcpp::IntegerVector location, int n, int threads) {
  const int N = X.n_rows, b = X.n_cols;
  if (location.size() != N) Rcpp::stop("Location map does not match the block");
  for (int i = 0; i < N; ++i)
    if (location[i] < 1 || location[i] > n) Rcpp::stop("Invalid location index");
  arma::mat S(n, b, arma::fill::zeros);
  arma::vec ss(b, arma::fill::zeros);
  const int *loc = location.begin();
#ifdef _OPENMP
  if (threads < 1) threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (int j = 0; j < b; ++j) {
    const double *x = X.colptr(j);
    double *s = S.colptr(j), acc = 0;
    for (int i = 0; i < N; ++i) {
      s[loc[i] - 1] += x[i];
      acc += x[i] * x[i];
    }
    ss[j] = acc;
  }
  return Rcpp::List::create(Rcpp::_["sums"] = S, Rcpp::_["ss"] = arma::accu(ss));
}
