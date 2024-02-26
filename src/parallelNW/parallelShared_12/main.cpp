#include <iostream>
#include <fstream>

#include <nwMpi.hpp>
#include <nwHelper.hpp>
#include <nwSolver.hpp>
#include "../../timer/timer.h"

#define TIME

using namespace NW_PAR;

int main(int argc, char const *argv[])
{
    MPI_Init(NULL, NULL);

#ifdef TIME
    if (argc < 7 || argc > 8) {
        std::cerr << "Program expects six or seven arguments:\nsequences (file), result (file), total time (file), \
                number of runs (integer), number of processes (integer), size of computation block (integer) [ time on process (file, optional) ]\n";
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    std::ofstream total_time(argv[3], std::fstream::app);
    const int runs = atoi(argv[4]);
    const int num_processes = atoi(argv[5]);
    const int block = atoi(argv[6]);
    std::ofstream process_time(argv[7], std::fstream::app);

    // declare variables for storing time
    myInt64 T_start_total, T_start_process;

    myDouble64 T_total, T_process;

    for(unsigned run = 0; run < runs + 4; run++) {
        if (argc == 8) T_process = 0;
#else
    if (argc != 3) {
        std::cerr << "Program expects two arguments: sequence, result (files)\n";
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
#endif


    MpiEnvironment env;
    NodeList nodes;
    NwMatrix matrix;
    seq_t short_seq;
    seq_t long_seq;

    prepareMpi(env, nodes);

    //std::cerr << "rank " << env.global_rank << " ready\n";

    if(env.global_rank == 0) {
        std::ifstream sequences_in;
        sequences_in.open(argv[1], std::ios::in);

        readSequences(sequences_in, short_seq, long_seq);

        sequences_in.close();

#ifdef TIME
        if (argc == 7) T_start_total = start_tsc();

#endif

    }

    shareSequences(env, nodes, matrix, short_seq, long_seq);

//    std::cerr << "rank " << env.global_rank << " got sequences\n";

    shareMatrix(env, matrix);

//    std::cerr << "rank " << env.global_rank << " got matrix\n";

#ifdef TIME
    env.block_size = block;
#else
    calculateBlockSize(env, matrix);
#endif

//    std::cerr << "rank " << env.global_rank << " got blocksize\n";

    setupScoreMatrix(matrix);

//    std::cerr << "rank " << env.global_rank << " prepared matrix\n";

    for (len_t start = 1; start < matrix.cols; start += env.block_size) {

        len_t end = std::min(start + env.block_size, matrix.cols);

        getPrevRow(env, matrix, start - 1, end);

#ifdef TIME
        if (argc == 8) T_start_process = start_tsc();
#endif

        computeBlock(matrix, short_seq, long_seq, start, end);

#ifdef TIME
        if (argc == 8) T_process += stop_tsc(T_start_process);
#endif

        sendLastRow(env, matrix, start - 1, end);

    }

//    std::cerr << "rank " << env.global_rank << " finished computation\n";

#ifdef TIME
    if (argc == 8 && run > 3) process_time << matrix.cols-1 << "\t" << T_process << "\t" << env.global_rank << std::endl;
#endif

    //syncronise al processes
    /*
     * int MPI_Barrier(MPI_Comm comm)
     */
    MPI_Barrier(env.global);

#ifdef TIME
   if (argc == 7 && env.global_rank == 0) {
        T_total = stop_tsc(T_start_total);
        if (run > 3) total_time << matrix.cols-1 << "\t" << T_total << "\t" << num_processes << std::endl;
    }
#endif

    if (!hasSuccessorOnNode(env) && !hasSuccessorOnOtherNode(env)) {
        score_t final_score = matrix.scores[matrix.rows*matrix.cols - 1];

        std::ofstream score_out;
        score_out.open(argv[2], std::ios::out);
        score_out << "score: " << final_score << std::endl;
        score_out.close();

    }


    cleanMpi(env, matrix, short_seq, long_seq);

#ifdef TIME
    } // end for (runs)
#endif

    MPI_Finalize();
    return 0;
}

