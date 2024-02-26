#include <userFunctions.hpp>
// functions to calculate gap penalties and similarity between bases

namespace NW_SEQ {
    score_t gapPenalty(const len_t gapLength) {
        // each gap gets the same penalty, regardless of the length
        return 0;

        // Small gaps are expensive, but do not increase much in
        // price as their length increases
        //if(gapLength == 0) return -4;
        //return -(3/gapLength); 
        
    }
    
    score_t similarity(const base_t a, const base_t b) {
        // positive score for match, negative for mismatch
        if ( a == b) return 1;
        else return -1;
    }

} // namespace NW_SEQ