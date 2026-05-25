#!/bin/bash
# Download reference datasets from CellxGene.
# These are used in steps 02-04 (reference preprocessing and label transfer).

wget --output-document all_neurons.h5ad \
    https://datasets.cellxgene.cziscience.com/0bb62fec-0cf1-46e1-9d10-de65a6d4f814.h5ad

wget --output-document all_non_neuronal.h5ad \
    https://datasets.cellxgene.cziscience.com/9343817c-5c74-4548-a7d5-01990df5af7a.h5ad
