#! /bin/bash

NTHREADS=$4
EXEC_ROOT=/lustre/scratch/zz374/CAMI/NB.run
GEN_SRC=$1
OUTPUT=$2
KMERSIZE=$3

MEMLIM=$(($NTHREADS * 1500000)) # 3 GB/core - 1.5 GB/class

$EXEC_ROOT train $GEN_SRC/ -s $OUTPUT -t $NTHREADS -k $KMERSIZE -m $MEMLIM
