#ifndef MAPT_INONETRACE_
#define MAPT_INONETRACE_
#include <madata.h>
#include <iostream>
namespace masb{
	struct pt_in_oneTrace_pram{
		std::string method = "fixe radius";
        int k_neighbours;
        float searchRadius;
    };
    void pt_in_oneTrace(pt_in_oneTrace_pram &pramas,PointList &ma_coords, VectorList &bisectors);
}

#endif