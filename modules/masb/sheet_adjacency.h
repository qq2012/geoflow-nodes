#ifndef SHEET_ADJACENCY_
#define SHEET_ADJACENCY_

#include <iostream>
#include "madata.h"

namespace ridge {
    struct adjacency_parameters {
        float searchRadius ;
        int  adjacency_thresh;
    };
    void adjacencyProcessing(ridge::adjacency_parameters& power, masb::MAT &mat, 
        masb::Sheet_idx_List &sheets, int_pair_vec &adjacency, masb::PointList &junction);
}
#endif // !SHEET_ADJACENCY_


