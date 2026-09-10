# SPDE timing

## README convex-boundary run (2026-09-10)

The US example was refitted on its existing synthetic BED with the new
coarse convex country envelope, buffered by a fixed 50 km independently of range.
The range is 290.7879 km and maximum triangle area
8,000 km². There are 1,213 mesh vertices, 872 supported vertices and rank 861;
`auto` selected `moments`. Nominal correlation is exp(-1) at 150 km.

| Stage | New SPDE seconds | Completed GP reference seconds |
| --- | ---: | ---: |
| Preparation | 2.03 | 3.49 |
| BED accumulation | 24.64 | 198.78 |
| REML | <0.01 | 0.01 |
| 10 PCs | 0.64 | 0.39 |
| Analysis total | **27.31** | **202.67** |

The reference GP run below was reused, not rerun. The regenerated BED, BIM,
FAM and frequency files have exactly the same MD5 hashes as its recorded inputs.
Both fit all 150,000 SNPs at 3,000 locations for 8,953 individuals, with the same
20,136,705 missing calls and fixed-effect columns. BED MD5 is
`5bc483bca91188fec4f0061f7c004ca5`; frequency-file MD5 is
`94e024650467282221d3374dd8117d55`. Generation and plotting are excluded.

The first three PCs recover the planted spatial subspace with canonical
correlations 0.98438, 0.97370, 0.93769 (new SPDE), compared with
0.98398, 0.97167, 0.93065 (GP reference). These permit sign changes and rotations.
Total time is 86.5% lower for SPDE in this comparison. This is a single rerun
against a completed reference, not a repeated or simultaneous benchmark.
GP rho is 135.853 km; this SPDE example uses 150 km as its e-folding distance.
The distinct kernels and mesh settings preclude a general accuracy ranking.

## Earlier matched-scale timing comparison

One sequential comparison on Windows 11, R 4.6.1, spatialGPCA 0.3.0, with
CppMatrix available and the default matrix-product backend. Both methods used
the same synthetic BED, external allele frequencies, 256-SNP blocks, 8,953
matched individuals, 3,000 locations, 150,000 SNPs and five fixed-effect columns
including the intercept. Both processed 20,136,705 missing calls and returned
10 orthonormal raw-location PCs.

| Stage | SPDE seconds | GP seconds |
| --- | ---: | ---: |
| Preparation | 2.13 | 3.49 |
| BED accumulation | 27.36 | 198.78 |
| REML | 0.02 | 0.01 |
| 10 PCs | 0.67 | 0.39 |
| Analysis total | **30.18** | **202.67** |

SPDE reduced accumulation time by 86.2% and analysis time by 85.1%. Mesh setup
and metadata added 0.11 seconds. Data generation, plotting and rendering are
excluded for both methods. Sparse B'WG took 1.54 seconds; SPDE's remaining
accumulation time was mainly covariance updates (15.29 s) and decoding (7.15 s).

The GP had 864 centres and rank 859. An independent RTriangle mesh had 1,416
vertices, 865 supported by observations and rank 856. Its maximum triangle area
was 8068.542 km² and buffer 135.853 km. The nominal Matérn correlation matched
exp(-1) at the GP rho of 135.853 km. This matches one physical scale, not the
full covariance model. SPDE selected the default `moments` accumulation path.
The US walkthrough uses rounded mesh/scale values; its images demonstrate the
workflow and are not this exact timing run.

For a larger mesh (7,087 total vertices, rank 4,276) with 15,000 locations and
2,000 SNPs, SPDE used 316.06 s for preparation, 27.49 s for accumulation and
137.85 s for 20 PCs: 481.41 s total including REML. Peak working set was 3.57 GiB.
This used the more stable `directions` accumulation path. Large-mesh preparation
and PCA remain dense costs; sparse projection does not guarantee lower total
time or memory for every dataset.
