/** 
 *  @file   nwSolver.hpp 
 *  @brief  Finds optimal alignment
 *  @author Sequence Aligners 
 *  @date   2021-11-15 
 ***********************************************/

#pragma once

#include <string>
#include <iostream>
#include <fstream>

#include <nwHelper.hpp>

namespace NW_SEQ_OPT {
    /**
     *  @brief Optimal Alignment NW
     *
     *  @details
     *   Finds an optimal alignment of two DNA sequences via Needleman Wunsch
     *   Fills score matrix with NW algorithm, does backtracking and stores outputs.
     *  @see   https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
     *  @todo  Implement correctness test with lots of test cases.
     * 
     *  @param inStream  to retrieve sequences
     *  @param outStream to write results to file
     *  @return     nothing
     */
    void findOptimalAlignment(NWMatrix& nw_matrix,const seq_t& seq1, const seq_t& seq2);

    
    
} // namespace NW_SEQ_OPT
