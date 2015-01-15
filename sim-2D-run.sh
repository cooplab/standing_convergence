#!/bin/bash

# Simulate time to mutation in a new patch in 2D
SB=".1"
N="100"
MU="1e-6"
THISDIR="sims"
mkdir -p $THISDIR
Rscript generate-patchy-run.R \
    "outdir=\"${THISDIR}\"" \
    "N=${N}" "range=c(101,101)" "sb=${SB}" "nsteps=40" "stepsize=50" "ntypes=10" "mu=${MU}"\
    's=matrix(sb,nrow=range[1],ncol=range[2])' \
    'initdist=array(c(rep(N,prod(range)),rep(0,(ntypes-1)*prod(range))),c(range,ntypes))' 
