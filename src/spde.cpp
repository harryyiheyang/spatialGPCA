#include "native.h"
#include <fstream>

// Triangle bounding boxes indexed by a uniform grid; no sf/fmesher runtime.
// [[Rcpp::export]]
Rcpp::List cpp_spde_project(const arma::mat &xy, const arma::imat &tv,
                          const arma::mat &loc) {
  if (xy.n_cols != 2 || tv.n_cols != 3 || loc.n_cols != 2 || !tv.n_rows ||
      !xy.is_finite() || !loc.is_finite() || tv.min() < 1 || tv.max() > (int)xy.n_rows)
    Rcpp::stop("Invalid planar mesh or locations");
  arma::rowvec lo = arma::min(xy, 0), span = arma::max(xy, 0) - lo;
  if (span.min() <= 0) Rcpp::stop("Mesh has zero extent");
  int ng = std::max(1, (int)std::sqrt(tv.n_rows));
  auto bin = [&](double z, int d) {
    return std::max(0, std::min(ng - 1, (int)std::floor((z - lo[d]) / span[d] * ng)));
  };
  std::vector<std::vector<int>> bins(ng * ng);
  const double tol = 1e-10;
  for (uword k = 0; k < tv.n_rows; ++k) {
    arma::mat z(3, 2);
    for (int j = 0; j < 3; ++j) z.row(j) = xy.row(tv(k,j) - 1);
    double det = (z(1,0)-z(0,0))*(z(2,1)-z(0,1)) -
                 (z(1,1)-z(0,1))*(z(2,0)-z(0,0));
    if (!std::isfinite(det) || det == 0) Rcpp::stop("Degenerate mesh triangle");
    for (int a = bin(z.col(0).min()-tol*span[0],0); a <= bin(z.col(0).max()+tol*span[0],0); ++a)
      for (int b = bin(z.col(1).min()-tol*span[1],1); b <= bin(z.col(1).max()+tol*span[1],1); ++b)
        bins[a + ng*b].push_back(k);
  }
  Rcpp::IntegerVector ii(3*loc.n_rows), jj(3*loc.n_rows);
  Rcpp::NumericVector ww(3*loc.n_rows);
  for (uword i = 0; i < loc.n_rows; ++i) {
    if (i % 4096 == 0) Rcpp::checkUserInterrupt();
    bool found = false;
    for (int k : bins[bin(loc(i,0),0) + ng*bin(loc(i,1),1)]) {
      int a = tv(k,0)-1, b = tv(k,1)-1, c = tv(k,2)-1;
      double ux=xy(b,0)-xy(a,0), uy=xy(b,1)-xy(a,1);
      double vx=xy(c,0)-xy(a,0), vy=xy(c,1)-xy(a,1);
      double dx=loc(i,0)-xy(a,0), dy=loc(i,1)-xy(a,1), det=ux*vy-uy*vx;
      double w[3]; w[1]=(dx*vy-dy*vx)/det; w[2]=(ux*dy-uy*dx)/det; w[0]=1-w[1]-w[2];
      if (std::min({w[0],w[1],w[2]}) < -tol) continue;
      double sum=0; for (double &v : w) { v=std::max(0.0,v); sum+=v; }
      for (int j=0;j<3;++j) { ii[3*i+j]=i+1; jj[3*i+j]=tv(k,j); ww[3*i+j]=w[j]/sum; }
      found=true; break;
    }
    if (!found) Rcpp::stop("A location is outside the SPDE mesh");
  }
  return Rcpp::List::create(Rcpp::_["i"]=ii,Rcpp::_["j"]=jj,Rcpp::_["x"]=ww);
}

struct SpdeStats {
  arma::sp_mat B;
  mat Xq, H, cov, projection;
  vec weights, rootw;
  double energy=0;
  uint64_t p=0, missing=0;
  std::map<std::string,double> timing;
  SpdeStats(const arma::sp_mat &b, const arma::mat &x, const arma::vec &w, const arma::mat &v) :
    B(b), Xq(x), projection(v), weights(w), rootw(arma::sqrt(w)) {
    if (!B.n_cols || B.n_rows != w.n_elem || x.n_rows != w.n_elem ||
        !w.is_finite() || arma::any(w<=0) || !x.is_finite())
      Rcpp::stop("Invalid SPDE accumulation dimensions or counts");
    H=B.t()*(Xq.each_col()%rootw);
    if (projection.n_cols && projection.n_rows != B.n_cols)
      Rcpp::stop("Invalid SPDE direction projection");
    uword dim = projection.n_cols ? projection.n_cols : B.n_cols;
    cov.zeros(dim,dim);
  }
  void add(mat G) {
    if (G.n_rows!=B.n_rows || !G.n_cols || !G.is_finite())
      Rcpp::stop("Invalid standardized location genotype block");
    auto t=Clock::now();
    mat U=B.t()*(G.each_col()%weights);
    timing["sparse_BtWG"]+=seconds(t);
    t=Clock::now();
    G.each_col()%=rootw;
    mat D=Xq.t()*G;
    U-=H*D;
    G-=Xq*D;
    energy+=arma::accu(arma::square(G));
    timing["fixed_effects_energy"]+=seconds(t);
    if (projection.n_cols) {
      t=Clock::now();
      U=projection.t()*U;
      timing["stable_direction_transform"]+=seconds(t);
    }
    t=Clock::now();
    int n=U.n_rows,k=U.n_cols; double one=1; char lower='L',nt='N';
    svgpc_dsyrk(&lower,&nt,&n,&k,&one,U.memptr(),&n,&one,cov.memptr(),&n);
    timing[projection.n_cols ? "direction_covariance" : "raw_covariance"]+=seconds(t);
    p+=G.n_cols;
  }
  Rcpp::List result() {
    return Rcpp::List::create(Rcpp::_["cross_covariance"]=arma::symmatl(cov),
      Rcpp::_["energy"]=energy,Rcpp::_["p"]=(double)p,
      Rcpp::_["missing"]=(double)missing,Rcpp::_["timing"]=timing,
      Rcpp::_["covariance_space"]=projection.n_cols ? "directions" : "moments");
  }
};

// [[Rcpp::export]]
Rcpp::List cpp_spde_accumulate_matrix(const arma::sp_mat &B, const arma::mat &Xq,
  const arma::vec &counts, const arma::mat &G, int block_size, const arma::mat &projection) {
  if (block_size<1 || !G.n_cols) Rcpp::stop("Empty SNP matrix or invalid block size");
  SpdeStats s(B,Xq,counts,projection);
  for (uword j=0;j<G.n_cols;j+=block_size) {
    Rcpp::checkUserInterrupt();
    s.add(G.cols(j,std::min(G.n_cols-1,j+block_size-1)));
  }
  return s.result();
}

// Reader contract: consecutive 1-based inclusive SNP indices -> n x block matrix.
// [[Rcpp::export]]
Rcpp::List cpp_spde_accumulate_reader(const arma::sp_mat &B, const arma::mat &Xq,
  const arma::vec &counts, Rcpp::Function reader, int p, int block_size, const arma::mat &projection) {
  if (p<1 || block_size<1) Rcpp::stop("Invalid SNP count or block size");
  SpdeStats s(B,Xq,counts,projection);
  for (int j=0;j<p;) {
    Rcpp::checkUserInterrupt();
    int end=j+std::min(block_size,p-j);
    auto t=Clock::now();
    mat G=Rcpp::as<mat>(reader(j+1,end));
    if (G.n_cols!=(uword)(end-j)) Rcpp::stop("Reader returned the wrong number of SNPs");
    s.timing["read_blocks"]+=seconds(t);
    s.add(std::move(G)); j=end;
  }
  return s.result();
}

// Same BED allele coding, missing imputation and location averaging as accumulate.cpp.
// [[Rcpp::export]]
Rcpp::List cpp_spde_accumulate_bed(const arma::sp_mat &B, const arma::mat &Xq,
  const arma::vec &counts, std::string path, Rcpp::IntegerVector locations,
  const arma::vec &af, int block_size, const arma::mat &projection) {
  if (block_size<1 || !af.n_elem || !af.is_finite() || arma::any(af<=0) || arma::any(af>=1))
    Rcpp::stop("Require positive block size and frequencies strictly between zero and one");
  SpdeStats s(B,Xq,counts,projection);
  vec check(counts.n_elem,arma::fill::zeros);
  for (int v : locations) {
    if (v<0 || v>(int)counts.n_elem) Rcpp::stop("Invalid sample to location mapping");
    if (v) check[v-1]++;
  }
  if (!arma::approx_equal(check,counts,"absdiff",0)) Rcpp::stop("Location counts differ from matched FAM samples");
  std::ifstream file(path,std::ios::binary);
  unsigned char header[3]; file.read(reinterpret_cast<char*>(header),3);
  if (!file || header[0]!=0x6c || header[1]!=0x1b || header[2]!=1)
    Rcpp::stop("Require a SNP-major PLINK BED file");
  size_t ns=locations.size(),stride=(ns+3)/4;
  file.seekg(0,std::ios::end);
  if (file.tellg()!=static_cast<std::streamoff>(3+stride*af.n_elem))
    Rcpp::stop("BED size disagrees with FAM or BIM");
  file.seekg(3,std::ios::beg);
  std::vector<unsigned char> bytes(stride);
  for (uword start=0;start<af.n_elem;start+=block_size) {
    Rcpp::checkUserInterrupt();
    uword bs=std::min((uword)block_size,af.n_elem-start);
    mat G(counts.n_elem,bs,arma::fill::zeros);
    auto t=Clock::now();
    for (uword j=0;j<bs;++j) {
      uword v=start+j; file.read(reinterpret_cast<char*>(bytes.data()),stride);
      if (!file) Rcpp::stop("Truncated BED");
      for (size_t i=0;i<ns;++i) {
        if (!locations[i]) continue;
        unsigned code=(bytes[i/4]>>(2*(i%4)))&3;
        double x=code==0 ? 2.0 : (code==2 ? 1.0 : 0.0);
        if (code==1) { s.missing++; x=2*af[v]; }
        G(locations[i]-1,j)+=x;
      }
      G.col(j)=(G.col(j)/counts-2*af[v])/std::sqrt(2*af[v]*(1-af[v]));
    }
    s.timing["decode_aggregate_standardize"]+=seconds(t);
    s.add(std::move(G));
  }
  return s.result();
}
