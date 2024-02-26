#!bin/bash
# This file takes on number N as input and generates sequences of length 2^i, where i goes from 1 to N.
# If no input is given, N = 10

#Set Output directories
OUTPUT_GEN=../data/sequences
#Remove old data
# Clean up

mkdir -p $OUTPUT_GEN #Create if it does not exist
rm -r $OUTPUT_GEN
mkdir $OUTPUT_GEN

#input defines number (and length of longest sequence) of generated sequences
N=$1
INIT_SIZE=${3:-100}
#second input defindes stepsize
STEP=${2:-10}

for n in $(seq $INIT_SIZE $STEP $N)
do
    #generate strings with increasing length
    #longest sequence has length n (other one can be shorter)
    python3 ../src/helperScripts/mutatedStrings.py $n > $OUTPUT_GEN/seq${n}.txt
done
