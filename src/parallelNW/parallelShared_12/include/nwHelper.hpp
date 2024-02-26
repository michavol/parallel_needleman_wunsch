#pragma once

/** 
 *  @file   nwHelper.hpp 
 *  @brief  Optimisations: For score matrix not matrix of structs
 *          but struct of matrices 
 *  @author Sequence Aligners 
 *  @date   2021-11-15 
 ***********************************************/

#include <cstdint>
#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>

namespace NW_PAR {
    /**
     * @brief datatype for encoding sequence elements
     *
     * The C++ Standard guarantees numbering 0, 1, 2, 3, 4, 5
     */
    enum base_t : unsigned char {A, T, C, G, NumBases, GAP};

    /**
     * @brief datatype for encoding predecessor of a cell in score matrix
     */
    enum pre_t : unsigned char {UP, DIAG, LEFT, UNDEF};

    /**
     * @brief datatype for scores
     */
    using score_t = int;

    /**
     * @brief length used for total length of sequences, gap lengths
     * 
     */
    using len_t = unsigned int;


    /**
     * @struct PreCollection - not needed at the moment
     * @brief a container for four predecessors
     *
     * one instance of this class must be of the same size as the datatypes for score and length
     *
     *
    struct PreCollection {
            static constexpr unsigned int size = sizeof(score_t)/sizeof(pre_t);
            pre_t data[size];
    };
    */


    /**
     * @struct NwLine
     * @brief store a single line of a NwMatrix
     *
     * the memory used for the arrays must be feed separately
     *
     */
    struct NwLine {
            len_t size;
            score_t* scores;
            len_t* gaps;
            pre_t* pres;
    };


    /**
     * @struct NwMatrix
     * @brief store score matrix
     *
     * store predecessor, gap length and score in three different matrices
     * matrices are contiguous arrays of dimension rows * cols
     *
     * the memory used for the arrays must be feed separately
     *
     * global row index does not include the first row that is filled with gaps
     *
     */
    struct NwMatrix {
            len_t global_row_index;
            len_t rows;
            len_t cols;
            score_t* scores;
            len_t* gaps;
            pre_t* pres;
            NwLine* prev_row;
    };


    /**
     * @brief container for sequences of bases
     *
     * the memory used for the arrays must be feed separately
     */
    struct seq_t{
            len_t size;
            base_t* data;
    };


    // /////////////////////////////////////////////////////////////////

    // INPUT FUNCTIONS - DEFINED IN SEPARATE FILE USERFUNCTIONS.CPP

    /**
     * @brief calculate gap penalty for a given gap length
     * @param[in] gap_length length of the gap
     * @return score of a gap of given length
     */
    score_t gapPenalty(len_t);

    /**
     * @brief compute score of a pair of bases that should be matched
     * @param[in] first base
     * @param[in] second base
     *
     * @return score of the pair
     */
    score_t similarity(base_t, base_t);

    // /////////////////////////////////////////////////////////////////

    /**
     * @brief transform a DNA sequence stored as a string to seq_t
     *
     * @param[in] s a DNA sequence
     * @param[out] seq the DNA sequence stored as seq_t
     */
    void stringToSequence(const std::string& s, seq_t& seq);


    /**
     * @brief read two DNA sequences from input
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
    void readSequences(std::istream& inStream,
            seq_t& seq1,
            seq_t& seq2);

    /**
     * @brief fill first column and row of global matrix with gaps
     * 
     * @param score_matrix memory for matrices is allocated in correct size
     * 
     * fill first column with gaps
     * if global row index is zero, fill prev_row
     *
     *
     */
    void setupScoreMatrix(NwMatrix& nw_matrix);

    /**
     * @brief
     *
     * @pre
     * @post
     * @param nw_matrix
     * @param seq1
     * @param seq2
     * @param start index of first column to be computed, must be strictly greater than 0
     * @param end past the end index of block to be computed, must be smaller or equal to nw_matrix.cols
     *
     * @todo implement
     */
    void computeBlock(NwMatrix& nw_matrix,
            const seq_t& seq1,
            const seq_t& seq2,
            const len_t start,
            const len_t end);

    /**
     * @brief backtrack & store optimal alignment
     * 
     * 
     * @param nw_matrix complete, filled score matrix
     * @param seq1 1srt sequence that should ne aligned
     * @param seq2 the sequences that are being aligned
     * 
     * @param aligned_seq1 reversed aligned sequences, including gaps
     * @param aligned_seq2 both sequences are of same size
     *
     * @return score_t total score of the alignment
     *
     * @todo implement
     */
    score_t backtrackSolution(const NwMatrix& nw_matrix,
            const seq_t& seq1, const seq_t& seq2,
            seq_t& aligned_seq1,
            seq_t& aligned_seq2);

    /**
     * @brief write total score and optimal alignment to outStream
     * 
     * @param outStream valid stream
     * @param total_score 
     * @param aligned_seq1 
     * @param aligned_seq2 
     *
     * @todo implement
     */
    void printAlignment(std::ostream& outStream,
            const score_t total_score,
            const seq_t& aligned_seq1,
            const seq_t& aligned_seq2);

    void printSequence(std::ostream& outStream, const seq_t& seq);

    void printScoreMatrix(std::ostream& outStream, const NwMatrix& nw_matrix);

    void printGapMatrix(std::ostream& outStream, const NwMatrix& nw_matrix);

    void printPreMatrix(std::ostream& outStream, const NwMatrix& nw_matrix);

    char convertToChar(base_t base);

    char convertToChar(pre_t pre);



} // namespace NW_PAR
