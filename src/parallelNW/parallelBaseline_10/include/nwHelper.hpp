#pragma once 

/**
 *  @file   nwHelper.hpp
 *  @brief  Helper Functions for general Sequence Alignment
 *  @author Sequence Aligners
 *  @date   2021-11-18
 ***********************************************/

#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>

namespace NW_PAR {

    // TYPE DEFINITIONS
    /**
     * @brief encoding for sequence elements
     *        The C++ Standard guarantees numbering 0, 1, 2, 3, 4, 5
     * 
     */
    enum base_t : unsigned char {A, T, C, G, NumBases, GAP};

    /**
     * @brief encoding predecessor of a cell in score matrix
     */
    enum pre_t : uint8_t {UP, DIAG, LEFT, UNDEF};

    /**
     * @brief datatype for scores
     */
    using score_t = int32_t;

    /**
     * @brief datatype for lengths and indices
     */
    using len_t = uint32_t;

    /**
     * @brief score matrix entries
     * 
     * @param score     of score 
     * @param gapLength must match the above defined types
     * @param predecessor 
     * 
     * changes in this datatype must be reflected in the definition of
     * the corresponding MPI_Datatype mpi_cell_t
     */ 
    struct Cell {
            int32_t score;
            uint32_t gapLength: 30;
            uint32_t predecessor: 2;
    };


    /**
     * @brief container for sequences of bases
     */
    using seq_t = std::vector<base_t>;

    /**
     * @brief container for score matrix
     */
    using score_mat_t = std::vector<std::vector<Cell> >;

    /**
     * @brief conainer for one line of scores
     */
    using score_line_t = std::vector<Cell>;



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

    // FUNCTIONS

    /**
     * @brief store two sequences from input
     * 
     * @param[in] inStream stream where sequences can be read
     * 
     * @param[out] seq1 first sequence, encoded according to base_t
     * @param[out[ seq2 second sequence, encoded according to base_t, is not shorter than seq1
     */
    void readSequences(std::istream& inStream,
            seq_t& seq1, seq_t& seq2);

    /**
     * @brief setup score matrix
     * 
     * @param[in, out] score_matrix container of size num_rows * num_cols, first column will be filled with gap penalty
     * 
     * @param[in] num_rows desired number of rows of score_matrix
     * @param[in] num_cols desired number of columns of score_matrix
     * @param[in] inital_gap_length  length of gap for first column, ie. number of rows above in full score matrix
     *
     * dimensions of score matrix may be greater than num_rows or num_cols
     */
    void initialiseScoreMatrix(score_mat_t& score_matrix,
            const len_t num_rows,
            const len_t num_cols,
            const len_t initial_gap_length);

    /**
     * @brief  fill a line with all gaps
     * @param[in, out] row to contain scores
     * @param[in] start_index first entry in row that should be filled with gap penalty
     * @param[in] num_entries number of entries that should be filled
     * @param[in] initial_gap_length length of the gap that is continued at position start_index
     */
    void fillRowWithGaps(score_line_t& row,
            const len_t start_index,
            const len_t num_entries,
            const len_t initial_gap_length);

    /**
     * @brief compute next cell
     *
     * if two cells have  same score, the following order is applied: UP - DIAGONAL - LEFT
     */
    Cell computeCell(const Cell& up, const Cell& diag, const Cell& left, base_t a, base_t b);

    /**
     * @brief compute part of the score matrix
     * 
     * @param[in] seq1 first sequence to align
     * @param[in] seq2 second sequence to align, at least as long as seq1
     *                  length of seq1 >= number of rows of score_matrix
     *                  length of seq2 == number of cols of score_matrix + 1
     * @param[in] row_above scores from the above row of the full score matrix
     * @param[in, out] score_matrix contains scores of alignment, filled up to column start_col - 1
     * @param[in] num_rows number of rows to be filled in score_matrix, starting from first row
     * @param[in] start_col column where computation of new scores should start
     * @param[in] num_cols number of columns that should be filled
     * 
     * use Needleman-Wunsch algorithm to compute the alignment scores of
     * full seq1 and seq2 in the interval [start_col, start_col + num_cols[
     *
     * size of score matrix may be bigger than length of sequences, computations start at the top row
     *
     **/
    void computeScoreBlock(const seq_t& seq1, const seq_t& seq2,
            const score_line_t& row_above,
            score_mat_t& score_matrix,
            const len_t num_rows,
            const len_t start_col,
            const len_t num_cols);


    /**
     * @brief backtrack and store optimal alignment
     * 
     * @param[in] score_matrix fully filled score matrix, first row is implicit
     * @param[in] seq1 first sequence that is to be aligned
     * @param[in] seq2 second sequence that is to be aligned
     * 
     * @param[out] aligned_seq1 aligned first sequence
     * @param[out] aligned_seq2 aligned second sequence
     * 
     * @return total score of the alignment
     *
     * @details aligned sequences are in reversed order. Both are of the same length and may include gaps
     *
     */ 
    score_t backtrackSolution(const score_mat_t& score_matrix,
            const seq_t& seq1,
            const seq_t& seq2,
            seq_t& aligned_seq1,
            seq_t& aligned_seq2);

    /**
     * @brief write total score and optimal alignment to outStream
     *
     * @param[in] outStream stream where alignment should be written
     * @param[in] total_score score of the alignment
     * @param[in] aligned_seq1 first sequence, aligned with second
     * @param[in] aligned_seq2 second sequence, aligned with first sequence
     *
     * @details aligned sequences must have the same length
     * the order of the bases is reversed compared to the unaligned sequences
     *
     */
    void printAlignment(std::ostream& outStream,
            const score_t total_score,
            const seq_t& aligned_seq1,
            const seq_t& aligned_seq2);
/**
 * @brief print score matrix 
 * 
 * @param score_matrix 
 * @return void 
 */
    void printScoreMatrix(const score_mat_t& score_matrix);

/**
 * @brief print Sequence 
 * 
 * @param sequence 
 */
    void printSequence(const seq_t& sequence);

} // namespace NW_PAR
