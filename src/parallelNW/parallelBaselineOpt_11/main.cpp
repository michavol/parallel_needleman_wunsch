#include <fstream>

#include <nwMpi.hpp>
#include "include/nwHelper.hpp"
#include "../../timer/timer.h"


#define TIME

using namespace NW_PAR;

 int main(int argc, char const *argv[])
 {
    // this process is doing i/o and backtracking
    const int head = 0;
    
    // setup MPI environment
    int my_rank;
    int comm_size;

    MPI_Init(NULL, NULL);

    // Create outstream
    #ifdef TIME
    std::ofstream out_benchmarking(argv[3], std::fstream::app);
    std::ofstream procs_occupation_benchmarking(argv[7], std::fstream::app);
    int runs = atoi(argv[4]);
    int procs = atoi(argv[5]);

    for(unsigned run = 0; run < runs + 4; run++){
    #endif

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    std::ifstream inStream;
    inStream.open(argv[1],std::ios::in);
    std::ofstream out_alignments(argv[2]);
    
    

    // perform algorithm

    // Get MPI information
    const int local_rank = my_rank;
    MPI_Comm communicator;
    MPI_Comm_dup(MPI_COMM_WORLD, &communicator);
    // decide how long a block should be for the computation
    #ifdef TIME
    const int block_length = atoi(argv[6]);
    #else
    const int block_length = 63;
    #endif
    

    // VARIABLE DECLARATION

    // sequences to be aligned
    seq_t short_sequence;
    seq_t long_sequence;

    // position in full score matrix, excluding the implicit first row
    // identically also the index of the first entry in the full short_sequence that is included in the local short_sequence
    len_t global_row_index;

    // dimensions of the local score matrix
    // num_rows does not include the first row, which is recieved from preceeding process
    // should be identical for all processes exept head
    len_t num_rows;

    // num_cols includes first column that is precomputed on initialisation
    len_t num_cols;

    // general height of score matrices
    // should be equal to num_rows on all processes exept head
    // this is relevant only for the head process to remember what size
    // score matrices of other processes have (actual height of head may be different)
    len_t block_height;

    // score matrix for each process
    //NWMatrix nw_matrix();

    // line that is recieved from preceeding process
    score_line_t score_row_above;
    gap_line_t gap_row_above;
    pred_line_t pred_row_above;
   

    // SETUP

    #ifdef TIME
    myInt64 T_start;
    #endif

    if (local_rank == head) {
        // read sequences from input
        readSequences(inStream, short_sequence, long_sequence);

        #ifdef TIME
        T_start = start_tsc();
        #endif

        // send sequences to other processes
        // num_rows and block_height are set in this function
        distributeSequences(short_sequence, long_sequence, num_rows, block_height,
                communicator, local_rank, comm_size);

        num_cols = long_sequence.size() + 1;
        global_row_index = 0;

        //initialise empty row
        score_row_above = score_line_t(num_cols);
        gap_row_above = gap_line_t(num_cols);
        pred_row_above = pred_line_t(num_cols);

        //initialise matrix for size of full score matrix
        //score_matrix = score_mat_t(short_sequence.size(), score_line_t(long_sequence.size() + 1));

        fillRowWithGaps(score_row_above, gap_row_above, pred_row_above, 0, num_cols, 0);

    } else {
        // get sequences from head
        // global_row_index is set in this functions

        recieveSequences(short_sequence, long_sequence, global_row_index,
                communicator, head, local_rank);

        num_rows = short_sequence.size();

        // length of sequence plus an initial row
        num_cols = long_sequence.size() + 1;

        block_height = num_rows;

        //initialise empty row
        score_row_above = score_line_t(num_cols);
        gap_row_above = gap_line_t(num_cols);
        pred_row_above = pred_line_t(num_cols);

        //initialise matrix  -> moved to shared scope
    }
    
    //initialise matrix with local sizes
    NWMatrix nw_matrix(num_rows, num_cols);

    //initialise score matrix, first column is precomputed as gaps
    //initial gap is equal to the global row index of this process plus one
    initialiseScoreMatrix(nw_matrix, num_rows, num_cols, global_row_index + 1);

    // COMPUTATION

    // iterate over blocks
    unsigned num_blocks = (num_cols - 1)/block_length;
    
    
    #ifdef TIME
    // Counter for total time spent in one processor
    myDouble64 cycles = 0.;
    myInt64 T_start_occupation;
    #endif

    for (unsigned i = 0; i < num_blocks; ++i) {
        //first column to be filled
        len_t current_col_index = i*block_length + 1;

        // setup row above
        // length of row above is always one larger than block_length
        // last entry of the block previously computed must be considered as well
        if (local_rank != head) {
            recieveSegmentOfRow(score_row_above, gap_row_above, pred_row_above, current_col_index - 1, block_length + 1, i,
                    communicator, local_rank-1, local_rank);
        }

        #ifdef TIME
        // Counter for total time spent in one processor
        if(argc == 8) T_start_occupation = start_tsc();
        #endif

        // compute block of score matrix
        computeScoreBlock(short_sequence, long_sequence,
                score_row_above, gap_row_above, pred_row_above, nw_matrix, num_rows, current_col_index, block_length);

        #ifdef TIME
        // Counter for total time spent on computation in one processor
        if(argc == 8) {
          myDouble64 T_stop = stop_tsc(T_start_occupation);
          cycles += T_stop;
        }
        #endif
    
        // send last row to successor
        if (local_rank < comm_size - 1) {
            sendSegmentOfRow(nw_matrix, num_rows - 1, current_col_index - 1,
                    block_length + 1, i, communicator, local_rank, local_rank + 1);

        }
 
    }
    

    // if length of longer sequence is not divisible by block_length
    len_t num_extra_cols = (num_cols - 1) % block_length;
    if (num_extra_cols > 0) {
        len_t current_col_index = num_blocks * block_length + 1;

        // setup row above
        if (local_rank != head) {
            recieveSegmentOfRow(score_row_above, gap_row_above, pred_row_above,current_col_index - 1, num_extra_cols + 1, num_blocks,
                    communicator, local_rank - 1, local_rank);
        }

        #ifdef TIME
        // Counter for total time spent in one processor
        if(argc == 8) T_start_occupation = start_tsc();
        #endif

        // compute block of score matrix
        computeScoreBlock(short_sequence, long_sequence,
                score_row_above, gap_row_above, pred_row_above, nw_matrix, num_rows, current_col_index, num_extra_cols);

        
        #ifdef TIME
        // Counter for total time spent in one processor
        if(argc == 8) {
          myDouble64 T_stop_occupation = stop_tsc(T_start_occupation);
          cycles += T_stop_occupation;
        }
        #endif

        // send last row to successor
        if (local_rank < comm_size - 1) {
            sendSegmentOfRow(nw_matrix, num_rows - 1, current_col_index - 1,
                    num_extra_cols + 1, num_blocks, communicator, local_rank, local_rank + 1);

        }
    }
    
    // syncronise
    MPI_Barrier(communicator);
    
    #ifdef TIME
    if (argc == 7 && local_rank == 0) {
         myDouble64 T_stop = stop_tsc(T_start);
         if(run > 3) out_benchmarking << num_cols-1 << "\t" << T_stop << "\t" << procs << '\t' << block_length << std::endl;
     }
    #endif
    if (local_rank == comm_size - 1) {
        score_t final_score = nw_matrix.scoreMatrix.back().back();
        out_alignments << "score = " << final_score << std::endl;
    }

    
    #ifdef TIME
    if(argc == 8 && run > 3) procs_occupation_benchmarking << long_sequence.size() << '\t' << cycles << '\t' << local_rank << std::endl;
    #endif
    
    // If we want processor number in runtimes data, then we need to oustream it
    
   
// /*
//     // BACKTRACKING - don't do bactracking
//     if (local_rank == head) {
//         // collect score matrix
//         recieveScoreMatrices(score_matrix, short_sequence.size(),
//                 long_sequence.size() + 1, num_rows, block_height,
//                 communicator, local_rank, comm_size);

//         #ifdef TIME
//         myDouble64 T_stop = stop_tsc(T_start);
//         if(run > 0) {
//             out_benchmarking << num_cols-1 << "\t" << T_stop << "\t" << procs << '\t' << block_length << std::endl;
//           }
//         #endif

//         // retrieve alignment
//         seq_t aligned_short;
//         seq_t aligned_long;
//         score_t score = backtrackSolution(score_matrix, short_sequence, long_sequence, aligned_short, aligned_long);

//         //print
//         printAlignment(out_alignments, score, aligned_short, aligned_long);
//     } else {
//         // send score matrix
//         sendScoreMatrix(score_matrix, num_rows, num_cols, global_row_index, communicator, local_rank, head);
//     }

//     // get in- and outstreams
//     // do alignment
//     //findOptimalAlignment(inStream, out_alignments, MPI_COMM_WORLD, my_rank, comm_size);
// */
    
    MPI_Comm_free(&communicator);
    #ifdef TIME
    }
    #endif

    // finalize
    MPI_Finalize();


   return 0;
}
 
