/** 
 *  @file   nwHelper.hpp 
 *  @brief  Helper Functions for Sequence Alignment
 *  @author Sequence Aligners 
 *  @date   2021-11-15 
 ***********************************************/

#pragma once 

#include <cstdint>
#include <iostream>
#include <vector>
#include <fstream>


namespace NW_SEQ_OPT {
    // sequence elements
    // Standard guarantees numbering 0, 1, 2, 3, 4, 5
    enum base_t : unsigned char {A, T, C, G, NumBases, GAP};
    
    // scores
    using score_t = int32_t;

    // length
    // used for total length of sequences, gap lengths
    using len_t = uint32_t;
        

    using matrix_score_t = std::vector<u_int32_t>;
    using matrix_gap_t = std::vector<u_int32_t>;
    using matrix_pred_t = std::vector<u_int8_t>; // Cannot be std::couted as integer!!

    using matrix_idx_t =  std::vector<std::vector<u_int32_t>>;

    /**
     * @brief score matrix entries
     * types of score and gapLength must match the above using
     * statements
     * changes in this datatype must be reflected in 
     * the definition of the corresponding MPI_Datatype mpi_cell_t
     */
    struct Cell {
        int32_t score;
        uint32_t gapLength: 30;
        uint32_t predecessor: 2; // 0: UP
                                 // 1 DIAG (UPPER LEFT)
                                 // 2 LEFT
                                 // 3 UNDEFINED
    };
    
    
    // container for sequences
    using seq_t = std::vector<base_t>;
    
    // container for score matrix
    using score_mat_t = std::vector<Cell>;
    
    // /////////////////////////////////////////////////////////////////

    /**
    * @brief store two DNA sequences from input
    * 
    * @param inStream a stream containing two sequences on separate lines
    *                 encoded as ASCII-Characters, must contain only characters 
    *                 used in base_t
    *                 
    * @param seq1 first sequence, encoded according to base_t
    * @param seq2 second sequence, encoded according to base_t
    * 
    * first sequence is always shorter than second sequence
    */
    void readSequences(std::ifstream& inStream,
                       seq_t& seq1, seq_t& seq2);

    /**
     * @brief setup score matrix
     * 
     * @param score_matrix    container for score matrix, dimensions must match
     * @param num_rows        define size of score matrix
     * @param num_cols        define size of score matrix
     * @param init_gap_length length of gap for first column, ie. number of
     *                        rows above this matrix in full score matrix 
     * 
     * @return matrix of size (num_rows + 1) * (num_cols + 1)
     *         first column and row are filled with gaps
     */
    void initialiseScoreMatrix(score_mat_t& score_matrix, const u_int64_t num_rows, const u_int64_t num_cols, const len_t init_gap_length);
                               
    /**
    * @brief: backtrack & store optimal alignment
     * 
     * 
     * @param score_matrix     complete, filled score matrix
     * @param seq1             1srt sequence that should ne aligned 
     * @param seq2             the sequences that are being aligned
     * 
     * @param aligned_seq1      reversed aligned sequences, including gaps
     * @param aligned_seq2      both sequences are of same size
     * @return score_t 
     */
    score_t backtrackSolution(const score_mat_t& score_matrix,
                              const seq_t& seq1, const seq_t& seq2,
                              seq_t& aligned_seq1,
                              seq_t& aligned_seq2);
    /**
     * @brief 
     * 
     * @param outStream write total score and optimal alignment to outStream
     * @param total_score 
     * @param aligned_seq1 
     * @param aligned_seq2 
     */
    
    void printAlignment(std::ofstream& outStream,
                        const score_t total_score,
                        const seq_t& aligned_seq1,
                        const seq_t& aligned_seq2);
    /**
     * @brief print score matrix with relevant information for debugging
     */
    void printScoreMatrix(const score_mat_t& score_matrix);
    
    
} // namespace NW_SEQ_OPT
