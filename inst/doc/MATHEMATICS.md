# QR-reduced spatial genetic PCA

Notation: n locations, p SNPs, m knots, W=diag(N_i), C=[1, X] with all
post-coordinate table columns averaged within matched locations. Let
U=orth(W^(1/2) C), q=rank(C), P=I-UU', Yw=W^(1/2)Y and Z=W^(1/2)B.
B_ij=exp(-||x_i-k_j||/rho), and Q is the
exponential knot kernel used as coefficient prior precision.

## Response regression and spatial QR reduction

Y1=P Yw is the weighted fixed-effect regression residual.
Compute a complete rank-revealing QR of A=Z'U. If s=rank(A), take
N=Q_A[,(s+1):m]. Then N has m-s columns, U'ZN=0, C'WBN=0.
The reduced design is ZN; this is NOT the m-column residualized design PZ.
After QR, applying P to ZN again is algebraically unnecessary and is omitted.

For each SNP:

    y = C beta + BN alpha + e
    alpha ~ N(0, sigma_gp^2 (N'QN)^(-1))
    e ~ N(0, sigma_e^2 W^(-1)), independent of alpha
    lambda = sigma_e^2 / sigma_gp^2

The user supplies rho. REML or GCV uses the genotype sufficient statistics
to select lambda.

## Fitting and PCA

Set QN=N'QN, S=(ZN)'ZN. Solve S V=QN V diag(d), V'QN V=I, retaining positive
eigenvalues. Then T=ZN V diag(d^(-1/2)), T'T=I, J=T'Y1 and

    alpha_hat = V diag(sqrt(d)/(d+lambda)) J
    theta_hat = N alpha_hat
    Phi = T / sqrt(W) = BN V diag(d^(-1/2))
    F_hat = BN alpha_hat = Phi diag(d/(d+lambda)) J

Division by sqrt(W) acts rowwise. The default raw_F_PCA returns orthonormal
left singular vectors of F_hat and squared singular values divided by p.
The H, H2 and 2H-H2 operators provide alternative spectral summaries.

## Noise independence and rank

Before aggregation, the working observational model has independent errors.
Location-mean errors have variances sigma_e^2/N_i, and W whitening gives
sigma_e^2 I_n. Regression residuals in n coordinates have covariance
sigma_e^2 P; they are correlated and singular. In any orthonormal n-q basis
K of U's orthogonal complement, K'e_w has covariance sigma_e^2 I_(n-q).
GP/noise independence is a zero cross-covariance assumption, not disjoint
column spaces or zero dot products between realized vectors.

Mean imputation and HWE scaling are fixed preprocessing rules. They do not
prove an exact iid model for real post-imputation genetic measurements; the
working count weights are retained, not replaced by per-SNP observed counts.

Let r=n-q, k=rank(ZN). Usually k=m-q when m>q and the designs have general
position. With m much smaller than n, k<r. With every distinct location used
as a knot and a nonsingular kernel, m=n and k=r can occur. This does not remove
observational noise or its independence from the spatial GP.

## Selection and the full-rank limit

E=||Y1||_F^2, a_i=||J_i||^2, b=E-sum(a_i), h_i=d_i/(d_i+lambda).
The relative profiled REML objective is

    sum(log(1+d/lambda)) + r log((b + sum(a*(1-h))) / E).

The relative log GCV objective is

    log((b + sum(a*(1-h)^2))/E) - 2 log((r-sum(h))/r).

These use a shared variance profile and pooled energy over SNPs.
No SNP-specific variance pair is returned.

When k=r, b is algebraically zero. Stable equivalent objectives, including
lambda=0, are

    REML = sum(log(d+lambda)) + r log(sum(a/(d+lambda))/E)
    GCV  = log(sum(a/(d+lambda)^2)/E)
           - 2 log(sum(1/(d+lambda))/r).

The powers of lambda cancel analytically. This avoids a spurious 0/0 or search
boundary error while retaining the same likelihood model. This is internal
optimizer support, not a request that users choose or fix lambda.
