#include "PtinDirectionLine.h"
#include <map>
using namespace ridge;
void PtinDirectionLine::processing(PtinDirectionLine_param& power, masb::PointList &pts, masb::intList& pt_id, masb::VectorList &vectors) {
    std::cout << "PtinDirectionLine::processing" << std::endl;

    std::map<int, int> seg_frequency;
    for (auto i : pt_id) {
        ++seg_frequency[i];
    }
    /*
    for (const auto& e : seg_frequency)
        std::cout << "Element " << e.first
         << " encountered " << e.second << " times\n";
    */
    for (const auto& e : seg_frequency) {
        auto cur_sheet = e.first;
        auto cur_size = e.second;
        if (cur_sheet == 0 || cur_size < 3)
            continue;

        std::cout << "cur_sheet ID = " << cur_sheet
            << "cur_size " << e.second << std::endl;

        masb::PointList cur_pt;
        cur_pt.reserve(cur_size);
        for (int i = 0; i < pt_id.size(); i++) {
            if (pt_id[i] == cur_sheet) {
                cur_pt.push_back(pts[i]);
            }
        }
        if (cur_pt.size() != cur_size)
            std::cout << "Error" << std::endl;


    }
}