using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Core
{
    class Transmissibility
    {
        public static double calculate(BaseBlock block_1, BaseBlock block_2)
        {
            // As each block stores a list of side faces areas, we need to find the index of each block relative to the other.
            int index_1, index_2;

            index_1 = Array.IndexOf(block_1.neighbour_blocks_indices, block_2.index);
            index_2 = Array.IndexOf(block_2.neighbour_blocks_indices, block_1.index);

            double length_1, area_1, permeability_1;
            double length_2, area_2, permeability_2;

            length_1 = block_1.delta_x_list[index_1];
            area_1 = block_1.area_list[index_1];
            permeability_1 = block_1.permeability_list[index_1];

            length_2 = block_2.delta_x_list[index_2];
            area_2 = block_2.area_list[index_2];
            permeability_2 = block_2.permeability_list[index_2];

            double G = 2 * Global.Bc / (length_1 / (area_1 * permeability_1) + length_2 / (area_2 * permeability_2));
            return G;
        }
    }
}
