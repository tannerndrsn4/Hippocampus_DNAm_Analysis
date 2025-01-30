#!/bin/bash

## This is the "submitR.sh" script used to run the PQLseq R script (runPQLseq.R) via a batch job 

ml R/4.4.2

Rscript Hippocampus_ARIMA_DNA.R
