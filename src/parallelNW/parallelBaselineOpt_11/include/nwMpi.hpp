#pragma once

/**
 *  @file   nwMpi.hpp
 *  @brief  Helper Functions to use MPI communication in Sequence Alignment
 *  @author Sequence Aligners
 *  @date   2021-11-18
 ***********************************************/

#include <mpi.h>

#include "nwHelper.hpp"

namespace NW_PAR {

    /**
     * Type declarations:
     * the following types are used in the send and recieve Functions
     * base_t <-> MPI_UNSIGNED_INT
     * Cell <-> MPI_LONG_LONG_INT
     *
     * changes in the desing for base_t and Cell must be reflected
     * in the types used for MPI communication
     */

    /**
     * @brief distribute sequences to all other ranks
     *
     * each process gets a segment of seq1
     * the remainder, i.e. the first local_block_height rows,
     * are treated on process local_rank
     * may use non-blocking sends, sequences should not be altered
     *
     * @param[in] seq1 first sequence, to be aligned with second sequence
     * @param[in] seq2 second sequence, not shorter than first sequence
     * @param[out] local_block_height length of the local part of first sequence
     * @param[out] other_block_height length of the part of first sequence for all other processes
     * @param[in] communicator MPI communicator that is used to exchange data
     * @param[in] local_rank rank of the process who calls the function
     * @param[in] comm_size number of processes in the communicator
     */
    void distributeSequences(const seq_t& seq1, const seq_t& seq2,
            len_t& local_block_height,
            len_t& other_block_height,
            MPI_Comm communicator,
            int local_rank, int comm_size);

    /**
     * @brief recieve sequences to align
     *
     * @param[out] seq1 first sequence, to be aligned with second sequence
     * @param[out] seq2 second sequence, not shorter than first sequence
     * @param[out] global_row_index index of first row of local score matrix in full score matrix, excluding the implicit first row
     * @param[in] communicator MPI communicator that is used to exchange data
     * @param[in] sender_rank rank of the process who sends the sequences
     * @param[in] local_rank rank of the process who calls the function
     */
    void recieveSequences(seq_t& seq1, seq_t& seq2,
            len_t& global_row_index,
            MPI_Comm communicator,
            int sender_rank, int local_rank);

    /**
     * @brief send part of generic Segment of row to the next process
     * 
     * @param nw_matrix    score_matrix local score matrix, columns up to end_col must be computed
     * @param row_index    start_col first entry that should be send
     * @param start_index  num_cols number of entries to be sent
     * @param num_entries  number of entries sent to score matrix
     * @param tag          tag identifies the message
     * @param communicator communicator MPI communicater that is used to exchange data
     * @param local_rank   local_rank rank of the process who calls the function
     * @param recieve_rank 
     */
    void sendSegmentOfRow(const NWMatrix& nw_matrix,
            const len_t row_index,
            const len_t start_index,
            const len_t num_entries,
            unsigned int tag,
            MPI_Comm communicator,
            int local_rank, int recieve_rank);


    /**
     * @brief recieve a line of score data from other process
     *
     * @param[out] row_above container containing the recieved scores
     * @param[in] start_index first position where data should be inserted in row_above
     * @param[in] num_entries number of data sets that should be recieved
     * @param[in] tag identifies the message
     * @param[in] communicator MPI communicator that is used to exchange data
     * @param[in] sender_rank rank of the process who sends the sequences
     * @param[in] local_rank rank of the process who calls the function
     */
    void recieveSegmentOfRow(score_line_t& score_row_above,
               gap_line_t& gap_row_above,
               pred_line_t& pred_row_above,
               const len_t start_index,
               const len_t num_entries,
               unsigned int tag,
               MPI_Comm communicator,
               int sender_rank, int local_rank);


    /**
     * @brief send score matrix to other process
     *
     * @param[in] score_matrix matrix to be sent
     * @param[in] num_rows number of rows in score_matrix
     * @param[in] num_cols number of columns in score_matrix
     * @param[in] communicator MPI communicator that is used to exchange data
     * @param[in] local_rank rank of the process who calls the function
     * @param[in] comm_size number of processes in the communicator
     */
    void sendScoreMatrix(const NWMatrix& nw_matrix,
            const len_t num_rows,
            const len_t num_cols,
            const len_t global_row_index,
            MPI_Comm communicator,
            int local_rank, int recieve_rank);

    /**
     * @brief collect full score matrix from contributions of other processes
     *
     * @param[out] score_matrix matrix containing all collected scores
     * @param[in] num_rows total number of rows of full score_matrix, excluding implicit first row
     * @param[in] num_cols total number of columns of full score_matrix
     * @param[in] local_block_height number of rows of local score matrix
     * @param[in] other_block_height number of rows of score matrix from other processes
     * @param[in] communicator MPI communicator that is used to exchange data
     * @param[in] local_rank rank of the process who calls the function
     * @param[in] comm_size number of processes in communicator
     */ 
    void recieveScoreMatrices(NWMatrix& nw_matrix,
            const len_t total_rows,
            const len_t total_cols,
            const len_t local_block_height,
            const len_t other_block_height,
            MPI_Comm communicator,
            int local_rank, int comm_size);



} // namespace NW_PAR
