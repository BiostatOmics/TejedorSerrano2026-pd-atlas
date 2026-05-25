"""Global configuration for scanpy-based scripts."""

import scanpy as sc

RANDOM_SEED = 0
SCANPY_N_JOBS = -1
SCANPY_VERBOSITY = 4

sc.settings.n_jobs = SCANPY_N_JOBS
sc.settings.verbosity = SCANPY_VERBOSITY
