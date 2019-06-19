#ifndef SHEET_ADJACENCY_
#define SHEET_ADJACENCY_

#include <iostream>
#include "madata.h"

namespace ridge {
    struct adjacency_parameters {
        float searchRadius ;
        int  adjacency_thresh;
        bool OnlySurfaceAdjacent;
        float DeepSurf_thresh;
    };
    void adjacencyProcessing(ridge::adjacency_parameters& power, masb::MAT &mat, masb::PointList &pointCloud,
        masb::Sheet_idx_List &sheets, int_pair_vec &adjacency, masb::PointList &junction, 
        masb::intList &junctionID);
}
#endif // !SHEET_ADJACENCY_


