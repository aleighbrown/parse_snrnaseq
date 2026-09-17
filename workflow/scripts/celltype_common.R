# =============================================================================
# Shared by celltype_prep.R, celltype_doublets.R and celltype_cluster.R:
# argument parsing, logging, threading. Sourced (not run) by those scripts,
# never by a worker process -- it parses the command line on load.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# --- arguments ---------------------------------------------------------------
# Every script is called as  Rscript <script> --key=value ...
parse_args <- function(args) {
  bad <- !grepl("^--[^=]+=", args)
  if (any(bad)) stop("Arguments must look like --key=value, got: ", paste(args[bad], collapse = " "))
  setNames(as.list(sub("^--[^=]+=", "", args)), sub("^--([^=]+)=.*", "\\1", args))
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

need_args <- function(keys) {
  miss <- setdiff(keys, names(opt))
  if (length(miss)) stop("Missing arguments: ", paste(miss, collapse = ", "))
}
num  <- function(k) as.numeric(opt[[k]])
flag <- function(k) tolower(opt[[k]]) %in% c("true", "1", "yes")
csv  <- function(k) { v <- trimws(strsplit(opt[[k]], ",")[[1]]); v[nzchar(v)] }
num_or_null <- function(k) if (tolower(opt[[k]]) %in% c("", "none", "null", "na")) NULL else num(k)

# --- logging -----------------------------------------------------------------
# The GB figure is R's own heap after a gc. Use it, and the peak from
# log_peak(), to size config/cluster.yaml.
t0 <- Sys.time()
step <- function(...) {
  mem <- sum(gc(verbose = FALSE)[, 2]) / 1024
  message(sprintf("[%s | %6.1f min | %6.1f GB] ", format(Sys.time(), "%H:%M:%S"),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), mem), ...)
}
log_peak <- function() {
  peak <- sum(gc(verbose = FALSE)[, 6]) / 1024
  step(sprintf("done. Peak R heap %.1f GB (main process only; worker processes and fread buffers come on top)", peak))
}

# --- threading ---------------------------------------------------------------
# NO fork() anywhere in these scripts. fread() runs OpenMP threads, and GNU
# OpenMP is not fork-safe: a forked child that touches OpenMP (xgboost inside
# scDblFinder did) can wait forever on a lock no thread will release. The
# first test run sat in scDblFinder's MulticoreParam until h_rt killed it.
# Parallel work uses fresh R processes (BiocParallel::SnowParam) instead.
#
# Fresh processes do NOT inherit set_blas(), and a Snow worker cannot call
# functions from this file either: worker functions pin their own threads
# (see azimuth_worker in celltype_cluster.R). Without that, OpenBLAS/xgboost
# in each scDblFinder worker used ~4 cores on a 2-slot job. The rules also
# export OMP_NUM_THREADS / OPENBLAS_NUM_THREADS = the job's slots.
set_blas <- function(n) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(n)
}
# BLAS, OpenMP (xgboost) and data.table all capped at n threads = the job's slots.
pin_threads <- function(n) {
  Sys.setenv(OMP_NUM_THREADS = n, OPENBLAS_NUM_THREADS = n)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(n)
    RhpcBLASctl::omp_set_num_threads(n)
  }
  data.table::setDTthreads(n)
}

setup_threads <- function(n) {
  setDTthreads(n)
  set_blas(1)                     # raised only around PCA
  options(future.globals.maxSize = Inf)
}

snow_param <- function(n, seed) {
  if (n > 1) {
    BiocParallel::SnowParam(workers = n, type = "SOCK", RNGseed = seed, log = TRUE)
  } else {
    BiocParallel::SerialParam(RNGseed = seed, log = TRUE)
  }
}

# RunUMAP reads its thread count (uwot threads, not processes) from the future
# plan. Only wrap calls that do NOT dispatch futures: under a multicore plan
# those would fork. ScaleData does dispatch futures, so it runs unwrapped.
with_workers <- function(n, expr) {
  if (n > 1 && future::supportsMulticore()) {
    old <- future::plan(future::multicore, workers = n)
    on.exit(future::plan(old), add = TRUE)
  }
  expr
}

# --- per-sublibrary files ----------------------------------------------------
# Written by celltype_prep.R, read by the other two. Names are the sublibrary
# names from config/sublibraries.tsv, which is what Snakemake expands over.
counts_file  <- function(dir, sublib) file.path(dir, "counts", paste0(sublib, ".rds"))
doublet_file <- function(dir, sublib) file.path(dir, paste0(sublib, ".tsv"))
