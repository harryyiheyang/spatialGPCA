# Run from the package source directory after installing spatialGPCA.
# Rscript tools/build_rmd.R [input.Rmd] [output_directory]
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1] else "README.Rmd"
destination <- if (length(args) > 1L) args[2] else NULL
rmarkdown::render(input, output_dir = destination,
                  envir = new.env(parent = globalenv()), encoding = "UTF-8")
