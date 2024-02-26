/** 
 *  @file   userFunctions.hpp 
 *  @brief  Functions to calculate gap penalties and similarity between bases
 *  @author Sequence Aligners 
 *  @date   2021-11-15 
 ***********************************************/

#pragma once 

#include <nwHelper.hpp>
// functions to calculate gap penalties and similarity between bases

namespace NW_SEQ_OPT {
    /* 
     * brief:
     * return gap penalty for given gapLength
     * 
     * input parameters:
     * gapLength    length of gap in current cell of user-defined similarity matrix
     * 
     * result:
     * return gap penalty
    */ 
    inline score_t gapPenalty(const len_t gapLength) {
        // each gap gets the same penalty, regardless of the length
        return 0;

        // Small gaps are expensive, but do not increase much in
        // price as their length increases
        //if(gapLength == 0) return -4;
        //return -(3/gapLength); 
        
    }

    /* 
     * brief:
     * return entry in user-defined similarity matrix
     * 
     * input parameters:
     * a, b         two bases to compare
     * 
     * result:
     * return similarity value
    */ 
    inline score_t similarity(const base_t a, const base_t b) {
        // positive score for match, negative for mismatch
        if ( a == b) return 1;
        else return -1;
    }
} // namespace NW_SEQ_OPT