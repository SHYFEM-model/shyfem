#!/bin/bash

bsub < shympi_seq.sh
bsub < shympi_mpi2.sh
bsub < shympi_mpi4.sh
bsub < shympi_mpi18.sh
bsub < shympi_mpi36.sh
bsub < shympi_mpi72.sh
bsub < shympi_mpi108.sh
bsub < shympi_mpi216.sh
