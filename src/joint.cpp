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

// X'Y for a tall N x b block X and tall N x l matrix Y, OpenMP over the columns of X.
// [[Rcpp::export]]
arma::mat cpp_joint_crossprod(const arma::mat &X, const arma::mat &Y, int threads) {
  const arma::uword N = X.n_rows, b = X.n_cols, l = Y.n_cols;
  if (Y.n_rows != N) Rcpp::stop("Nonconformable crossprod");
  arma::mat Yt = Y.t();  // l x N, row i of Y contiguous
  arma::mat out(l, b);
#ifdef _OPENMP
  if (threads < 1) threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (arma::uword j = 0; j < b; ++j) {
    std::vector<double> acc(l, 0.0);
    const double *x = X.colptr(j), *yt = Yt.memptr();
    for (arma::uword i = 0; i < N; ++i) {
      const double xi = x[i];
      if (xi == 0) continue;
      const double *yr = yt + i * l;
      for (arma::uword c = 0; c < l; ++c) acc[c] += xi * yr[c];
    }
    for (arma::uword c = 0; c < l; ++c) out(c, j) = acc[c];
  }
  return out.t();  // b x l
}

// X T for a tall N x b block X and small b x l matrix T, OpenMP over row chunks.
// [[Rcpp::export]]
arma::mat cpp_joint_prod(const arma::mat &X, const arma::mat &T, int threads) {
  const arma::uword N = X.n_rows, b = X.n_cols, l = T.n_cols;
  if (T.n_rows != b) Rcpp::stop("Nonconformable product");
  arma::mat Tt = T.t();  // l x b: row j of T contiguous
  arma::mat Ot(l, N, arma::fill::zeros);
  const arma::uword chunk = 512, nchunk = (N + chunk - 1) / chunk;
#ifdef _OPENMP
  if (threads < 1) threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (arma::uword h = 0; h < nchunk; ++h) {
    const arma::uword i0 = h * chunk, i1 = std::min(N, i0 + chunk);
    double *o = Ot.memptr();
    for (arma::uword j = 0; j < b; ++j) {
      const double *x = X.colptr(j), *t = Tt.colptr(j);
      for (arma::uword i = i0; i < i1; ++i) {
        const double xi = x[i];
        double *orow = o + i * l;
        for (arma::uword c = 0; c < l; ++c) orow[c] += xi * t[c];
      }
    }
  }
  return Ot.t();
}

// Fused scan of one BED block straight from the 2-bit codes, without forming the N x b double
// matrix. With x_ij = lut_j[code_ij] (standardized, missing = 0):
//   sums (n x b): location sums Z'X;  ss: sum of x^2;  XtY (b x c): X'Y with Y given as Yt (c x N);
//   RtT (l x N, updated in place) += (X T)' with T = XtY[, q:(c-1)], the subspace-iteration product.
// Phase 1 runs over SNPs, phase 2 over chunks of individuals; both use OpenMP.
// [[Rcpp::export]]
Rcpp::List cpp_joint_bed_scan(std::string path, int file_samples, Rcpp::IntegerVector rows,
                              Rcpp::IntegerVector snp, const arma::vec &af,
                              Rcpp::IntegerVector location, int n, const arma::mat &Yt, int q,
                              Rcpp::NumericMatrix RtT, int threads) {
  // snp: 0-based variant indices of this block (contiguous runs are read in one call); af: one
  // frequency per entry of snp.
  const int count = snp.size();
  if (file_samples < 1 || count < 1 || (int)af.n_elem != count) Rcpp::stop("Invalid BED block request");
  const int N = rows.size(), c = Yt.n_rows, l = c - q;
  if ((int)Yt.n_cols != N || location.size() != N || q < 0 || q > c) Rcpp::stop("Nonconformable scan inputs");
  if (l > 0 && (RtT.nrow() != l || RtT.ncol() != N)) Rcpp::stop("Accumulator has the wrong shape");
  const size_t stride = ((size_t)file_samples + 3) / 4;
  std::ifstream file(path, std::ios::binary);
  if (!file) Rcpp::stop("Cannot open BED file");
  unsigned char header[3];
  file.read(reinterpret_cast<char *>(header), 3);
  if (!file || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("Require a SNP-major PLINK BED file");
  std::vector<unsigned char> bytes(stride * (size_t)count);
  for (int j = 0; j < count;) {
    if (snp[j] < 0) Rcpp::stop("Invalid variant index");
    int e = j + 1;
    while (e < count && snp[e] == snp[e - 1] + 1) ++e;
    file.seekg(3 + (std::streamoff)stride * snp[j], std::ios::beg);
    file.read(reinterpret_cast<char *>(bytes.data() + stride * (size_t)j), (std::streamsize)(stride * (size_t)(e - j)));
    if (!file) Rcpp::stop("Truncated BED");
    j = e;
  }
  std::vector<int> r(N), loc(N);
  for (int i = 0; i < N; ++i) {
    if (rows[i] < 1 || rows[i] > file_samples) Rcpp::stop("Invalid FAM row index");
    if (location[i] < 1 || location[i] > n) Rcpp::stop("Invalid location index");
    r[i] = rows[i] - 1; loc[i] = location[i] - 1;
  }
  // Byte offset and shift of every analysed individual, shared by both phases.
  std::vector<size_t> off(N); std::vector<unsigned char> sh(N);
  for (int i = 0; i < N; ++i) { off[i] = (size_t)(r[i] >> 2); sh[i] = (unsigned char)(2 * (r[i] & 3)); }
  std::vector<double> lut((size_t)4 * count);
  for (int j = 0; j < count; ++j) {
    const double f = af[j], mu = 2 * f, sc = 1 / std::sqrt(2 * f * (1 - f));
    // PLINK codes: 00 hom A1 (2 copies), 01 missing, 10 het, 11 hom A2.
    lut[4 * j] = (2 - mu) * sc; lut[4 * j + 1] = 0; lut[4 * j + 2] = (1 - mu) * sc; lut[4 * j + 3] = -mu * sc;
  }
  arma::mat S(n, count, arma::fill::zeros), XtY(count, c);
  arma::vec ssj(count);
  const double *yt = Yt.memptr();
  // Phase 1 is tiled: SNPs in groups of 16 per task, individuals in chunks of 512, so a chunk of Yt
  // stays in cache across the group instead of streaming all of Yt once per SNP.
  const int group = 4, ngroup = (count + group - 1) / group, rchunk = 512;
#ifdef _OPENMP
  if (threads < 1) threads = omp_get_max_threads();
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
  for (int g = 0; g < ngroup; ++g) {
    const int j0 = g * group, j1 = std::min(count, j0 + group), gj = j1 - j0;
    // Sums of Y rows by genotype class (adds only); X'Y = sum over classes of lut * class sum.
    std::vector<double> acc((size_t)gj * 3 * c, 0.0), cnt((size_t)gj * 3, 0.0);
    for (int i0 = 0; i0 < N; i0 += rchunk) {
      const int i1 = std::min(N, i0 + rchunk);
      for (int j = j0; j < j1; ++j) {
        const unsigned char *b = bytes.data() + stride * (size_t)j;
        const double *L = lut.data() + 4 * j;
        double *a = acc.data() + (size_t)(j - j0) * 3 * c, *nc = cnt.data() + (size_t)(j - j0) * 3;
        double *s = S.colptr(j);
        for (int i = i0; i < i1; ++i) {
          const int code = (b[off[i]] >> sh[i]) & 3;
          if (code == 1) continue;
          const int cls = code == 0 ? 0 : code - 1;  // 0: hom A1, 1: het, 2: hom A2
          const double *y = yt + (size_t)i * c;
          double *ac = a + (size_t)cls * c;
          for (int h = 0; h < c; ++h) ac[h] += y[h];
          s[loc[i]] += L[code];
          nc[cls] += 1;
        }
      }
    }
    for (int j = j0; j < j1; ++j) {
      const double *L = lut.data() + 4 * j, *a = acc.data() + (size_t)(j - j0) * 3 * c,
                   *nc = cnt.data() + (size_t)(j - j0) * 3;
      for (int h = 0; h < c; ++h) XtY(j, h) = L[0] * a[h] + L[2] * a[c + h] + L[3] * a[2 * c + h];
      ssj[j] = nc[0] * L[0] * L[0] + nc[1] * L[2] * L[2] + nc[2] * L[3] * L[3];
    }
  }
  if (l > 0) {
    // Per-SNP tables lut_j[code] * T_j (4 x l), then RtT[, i] += table_j[code_ij] over SNPs.
    std::vector<double> tab((size_t)count * 4 * l);
    for (int j = 0; j < count; ++j)
      for (int code = 0; code < 4; ++code)
        for (int h = 0; h < l; ++h) tab[((size_t)j * 4 + code) * l + h] = lut[4 * j + code] * XtY(j, q + h);
    double *out = RtT.begin();
    const int chunk = 256, nchunk = (N + chunk - 1) / chunk;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int ch = 0; ch < nchunk; ++ch) {
      const int i0 = ch * chunk, i1 = std::min(N, i0 + chunk);
      for (int j = 0; j < count; ++j) {
        const unsigned char *b = bytes.data() + stride * (size_t)j;
        const double *T = tab.data() + (size_t)j * 4 * l;
        for (int i = i0; i < i1; ++i) {
          const int code = (b[off[i]] >> sh[i]) & 3;
          if (code == 1) continue;
          const double *t = T + code * l;
          double *o = out + (size_t)i * l;
          for (int h = 0; h < l; ++h) o[h] += t[h];
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::_["sums"] = S, Rcpp::_["ss"] = arma::accu(ssj), Rcpp::_["XtY"] = XtY);
}
