# SPDE timing

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
