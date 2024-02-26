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
     *  @brief Cell Data type with three matrices to improve regularity in loops
     *
     *  @details
     *   ...
     *  @todo  Implement
     */
    struct NWMatrix {
        // Constructor to initialize matrices with given size - Reserve space such that memory will be contiguous
        NWMatrix(u_int64_t rows, u_int64_t cols) : m(rows), n(cols) {
            scoreMatrix = matrix_score_t(m*n);
            gapMatrix = matrix_gap_t(m*n);
            predMatrix = matrix_pred_t(m*n);
        }

        // Struct members
        matrix_score_t scoreMatrix;
        matrix_gap_t gapMatrix;
        matrix_pred_t predMatrix;

        u_int64_t m;
        u_int64_t n;
    };
    
    // container for sequences
    using seq_t = std::vector<base_t>;
    

    // /////////////////////////////////////////////////////////////////

    /* 
     * brief:
     * store two DNA sequences from input
     * 
     * input parameters:
     * inStream     a stream containing two sequences on separate lines
     *              encoded as ASCII-Characters, must contain only 
     *              characters used in base_t
     * 
     * result:
     * inStream     two lines have been read
     * seq1         first sequence, encoded according to base_t
     * seq2         second sequence, encoded according to base_t
     * first sequence is always shorter than second sequence
    */ 
    void readSequences(std::ifstream& inStream,
                       seq_t& seq1, seq_t& seq2);

    /*
     * brief:
     * setup score matrix
     * 
     * input:
     * score_matrix     empty container
     * num_rows
     * num_cols         define size of score matrix
     * init_gap_length  length of gap for first column, ie. number of 
     *                  rows above this matrix in full score matrix
     * 
     * result:
     * score_matrix     matrix of size (num_rows + 1) * (num_cols + 1)
     *                  first column and row are filled with gaps
     */
    void initialiseScoreMatrix(NWMatrix& nw_matrix, const u_int64_t num_rows, const u_int64_t num_cols, const len_t init_gap_length);


    /*
     * brief:
     * backtrack & store optimal alignment
     * 
     * input:
     * score_matrix     complete, filled score matrix
     * seq1
     * seq2             the sequences that are being aligned
     * 
     * 
     * result:
     * aligned_seq1
     * aligned_seq2     reversed aligned sequences, including gaps
     *                  both sequences are of same size
     * 
     * returns: 
     * total score of the alignment
     */ 
    score_t backtrackSolution(const NWMatrix& nw_matrix,
                              const seq_t& seq1, const seq_t& seq2,
                              seq_t& aligned_seq1,
                              seq_t& aligned_seq2);

    /*
     * brief:
     * write total score and optimal alignment to outStream
     */
    void printAlignment(std::ofstream& outStream,
                        const score_t total_score,
                        const seq_t& aligned_seq1,
                        const seq_t& aligned_seq2);


    /*
     * brief:
     * print score matrix with relevant information for debugging
     */
    void printScoreMatrix(const NWMatrix& nw_matrix);
    
    
} // namespace NW_SEQ_OPT
