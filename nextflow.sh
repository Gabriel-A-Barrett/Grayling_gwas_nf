#!/usr/bin/sh

nextflow run main.nf -entry GWAS_NF -with-singularity gwas-nf.sif -profile slurm

