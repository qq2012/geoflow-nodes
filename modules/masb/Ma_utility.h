#ifndef MA_UTILITY_
#define MA_UTILITY_

#include "madata.h"
namespace masb{
	struct Filter8R{
    public:
        float radius = 200;	 	
	}; 
    class idx_filter {
    public:
        void processing(ma_data &madata, ma_Geometry &maGeometry,intList &remaining_idx, 
            ma_data &remainingData, intList &remainingma_in_out, ma_Geometry &remainingGeometry);
        void processing(ma_data &madata, intList &remaining_idx,
            ma_data &remainingData, intList &remainingma_in_out);
    };
}

#endif