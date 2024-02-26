#include "nwHelper.hpp"

// functions to calculate gap penalties and similarity between bases

namespace NW_PAR {
    score_t gapPenalty(len_t gapLength) {
        // each gap gets the same penalty, regardless of the length
        return 0;
    }
    
    score_t similarity(base_t a, base_t b) {
        // positive score for match, negative for mismatch
        if ( a == b) return 1;
        else return -1;
    }

} // namespace NW
