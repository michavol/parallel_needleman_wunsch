# DPHPC project

Sequence alignment is salient for solving many fundamental problems in bioinformatics and can be achieved using the well-known dynamic programming Needleman--Wunsch (NW) algorithm. In this project, several approaches for optimising this algorithm are proposed and compared. This includes considering four distinct data types, antidiagonal and row-wise ways of traversing the dynamic programming matrix, increasing cache locality via index transformation and the application of SIMD instructions. For parallelising the NW algorithm, two different concepts were considered. Multi-core OpenMP can be employed by leveraging the antidiagonal way of traversing the NW matrix and multi-node MPI is implemented by splitting the work into bands of blocks where each band can be executed by one process. Besides using different datatypes, this project further expands on the idea by using shared memory. Benchmarks in the directory 'legacy' were performed on the Euler IV cluster equipped with Xeon Gold 6150 processors (2.7-3.7 GHz), provided by ETH Zurich.

To run this project:

```
mkdir build
cd build
cmake ..
cd ..
sh benchmark.sh N STEP RUNS PROCS BLOCK_LENGTH INIT_SIZE
``` 
- N defines length of the longest sequence
- STEP defines stepsize
- RUNS defines number of iterations
- PROCS defines number of maximal processors used for parallel implementations
- BLOCK_LENGTH defines maximal block length
- INIT_SIZE defines the initial sequence length


To visualise runtimes:

```
sh visualization.sh N
```
- N defines length of the longest sequence

sequential implementations:

    - ALG00 Cell Matrix
    - ALG01 Integer Matrices
    - ALG02 Diagonal Index Matrix
    - ALG03 Diagonal Integer Vectors
    - ALG04 Diagonal Integer Vectors SIMD
    - ALG05 Diagonal Integer Vectors Nested
    - ALG06 Diagonal Integer Vectors Pragma SIMD

parallel implementations:

    - ALG10 MPI Cell Matrix
    - ALG11 MPI Integer Matrices
    - ALG12 Shared MPI Integer Vectors
    - ALG13 OpenMP Cell Vector
    - ALG14 OpenMP Integer Vector

dependencies:

    - biopython
    - matplotlib
    - numpy
    - seaborn
    - cmake
    - openMPI
