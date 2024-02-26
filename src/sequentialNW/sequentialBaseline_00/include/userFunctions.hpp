/**
 * @file userFunctions.hpp
 * @author Sequence Aligners
 * @brief Functions to calculate gap penalties and similarity between bases
 * @version 0.1
 * @date 2021-11-21
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once 

#include <nwHelper.hpp>

namespace NW_SEQ {
   /**
    * @brief return gap penalty for given gapLength
    * 
    * @param gapLength length of gap in current cell of user-defined similarity matrix
    * @return return gap penalty
    */
    score_t gapPenalty(const len_t gapLength);

   /**
    * @brief return entry in user-defined similarity matrix
    * 
    * @param a one base to compare
    * @param b one base to compare
    * @return return similarity value
    */
    score_t similarity(const base_t a, const base_t b);
} // namespace NW_SEQ