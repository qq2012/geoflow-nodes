#ifndef MA_UTILITY_
#define MA_UTILITY_

#include "madata.h"
namespace masb{
	struct Filter8R{
    public:
        float radius;	 	
	}; 
    class ind_filter {
    public:
        void processing(ma_data &madata, ma_Geometry &maGeometry,intList &remaining_idx);
    };
}

#endif