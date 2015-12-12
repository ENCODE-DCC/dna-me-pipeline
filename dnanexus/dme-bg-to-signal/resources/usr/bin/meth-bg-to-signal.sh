#!/bin/bash -e

if [ $# -ne 2 ]; then
    echo "usage v1: meth-bg-to-signal.sh <bedGraph> <chrom.sizes>"
    echo "Converts bismark2bedGraph output to bigWig file.  Is independent of DX and ENCODE."
    exit -1; 
fi
bedgraph=$1 # uncompressed bedGraph output from bismark2bedGraph.
chrom_sizes=$2 # chrom.sizes file used by bedGraphToBigWig

target_root=${bedgraph%.bedGraph}
target_root=${target_root%.bg}

echo "-- Convert to signal bedGraph to bigWig..."
set -x
bedGraphToBigWig $bedgraph $chrom_sizes ${target_root}.bw
set +x

echo "-- The results..."
ls -l ${target_root}*

