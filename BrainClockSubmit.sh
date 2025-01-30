#!/bin/bash

## This is the "submitR.sh" script used to run the PQLseq R script (runPQLseq.R) via a batch job 

module load R/4.3.2

Rscript Clock_construction_hippocampus.R
