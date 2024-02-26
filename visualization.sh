#!/bin/bash
# This file creates plots to compare runtimes 

OUTPUT_PNG=benchmarking/plots

# Clean up
mkdir -p $OUTPUT_PNG 
rm -rf $OUTPUT_PNG
mkdir $OUTPUT_PNG

N_SEQ=${1:-1000}
N_SEQ_P=${2:-$N_SEQ}
N_SEQ_P_OCC=${3:-$N_SEQ}
N_SEQ_BLOCK=${4:-$N_SEQ}

cd benchmarking

rm -rf plots
mkdir plots

# Set Input/Output directories
RUNTIMES=runtimes

# Define algorithm directory names
ALG_00=sequentialBaseline_00
ALG_01=sequentialOpt_01
ALG_02=sequentialDiag_02
ALG_03=sequentialDiagOpt_03
ALG_04=sequentialDiagSimd_04
ALG_05=sequentialDiagOptNested_05
ALG_06=sequentialDiagSimdPragma_06
ALG_10=parallelBaseline_10
ALG_11=parallelBaselineOpt_11
ALG_12=parallelShared_12 #included in the list, but not in the python scripts besides seqlength_vs_cycles
ALG_13=parallelOpenMp_13
ALG_14=parallelOpenMpOpt_14


# Create seqlength vs. cycles plot
# Remember to run the following before
# 'sh seqlength_vs_cycles.sh <Max length of sequences> <Step size of sequence lengths> <Runs for averaging for benchmarking> <Number Processors>'
#python3 plot_seqlength_vs_cycles.py $RUNTIMES/$ALG_00/runtime_${ALG_00}_${N_SEQ}.txt $RUNTIMES/$ALG_01/runtime_${ALG_01}_${N_SEQ}.txt $RUNTIMES/$ALG_02/runtime_${ALG_02}_${N_SEQ}.txt $RUNTIMES/$ALG_03/runtime_${ALG_03}_${N_SEQ}.txt $RUNTIMES/$ALG_04/runtime_${ALG_04}_${N_SEQ}.txt $RUNTIMES/$ALG_05/runtime_${ALG_05}_${N_SEQ}.txt $RUNTIMES/$ALG_06/runtime_${ALG_06}_${N_SEQ}.txt $RUNTIMES/$ALG_10/runtime_${ALG_10}_${N_SEQ}.txt $RUNTIMES/$ALG_11/runtime_${ALG_11}_${N_SEQ}.txt $RUNTIMES/$ALG_12/runtime_${ALG_12}_${N_SEQ}.txt $RUNTIMES/$ALG_13/runtime_${ALG_13}_${N_SEQ}.txt $RUNTIMES/$ALG_14/runtime_${ALG_14}_${N_SEQ}.txt 

# Create procs vs. cycles plot
# Remember to run the following before
# 'sh seqlength_vs_cycles.sh <Max length of sequences> <Step size of sequence lengths> <Runs for averaging for benchmarking> <Number Processors>'
python3 plot_procs_vs_cycles.py     runtimes_procs_vs_cycles/$ALG_10/runtimes_procs_vs_cycles_${ALG_10}_${N_SEQ_P}.txt  runtimes_procs_vs_cycles/$ALG_11/runtimes_procs_vs_cycles_${ALG_11}_${N_SEQ_P}.txt runtimes_threads_vs_cycles/$ALG_13/runtimes_threads_vs_cycles_${ALG_13}_${N_SEQ_P}.txt runtimes_threads_vs_cycles/$ALG_14/runtimes_threads_vs_cycles_${ALG_14}_${N_SEQ_P}.txt

# Create processor and thread occupation plot
#python3 plot_procs_occupation.py    runtimes_procs_occupation/$ALG_10/runtimes_procs_occupation_${ALG_10}_${N_SEQ_P_OCC}.txt    #runtimes_threads_occupation/$ALG_13/runtimes_threads_occupation_${ALG_13}_${N_SEQ_P_OCC}.txt

# Create block_length vs cycles plot
#python3 plot_blockLength_vs_cycles.py    runtimes_blockLength_vs_cycles/$ALG_10/runtimes_blockLength_vs_cycles_${ALG_10}_${N_SEQ_BLOCK}.txt
#python3 plot_blockLength_vs_cycles.py    runtimes_blockLength_vs_cycles/$ALG_11/runtimes_blockLength_vs_cycles_${ALG_11}_${N_SEQ_BLOCK}.txt



