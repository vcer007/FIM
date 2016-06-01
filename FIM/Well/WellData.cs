using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Well
{
    class WellData
    {
        public static void setR_equivalent(BaseBlock block, BaseBlock[] grid)
        {
            double r_equivalnt = 0;

            double numerator = 0, denominator = 0;

            BaseBlock neighbour;
            // neighbour block internal index of "block".
            int index;

            for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
            {
                if (block.neighbour_blocks_indices[i] >= 0)
                {
                    neighbour = grid[block.neighbour_blocks_indices[i]];
                    index = Array.IndexOf(neighbour.neighbour_blocks_indices, block.index);

                    numerator += block.boundary_length_list[i] / (block.delta_x_list[i] + neighbour.delta_x_list[index]) * Math.Log(block.delta_x_list[i] + neighbour.delta_x_list[index]);
                    denominator += (block.boundary_length_list[i]) / (block.delta_x_list[i] + neighbour.delta_x_list[index]);
                }
            }

            numerator -= Global.PI;

            r_equivalnt = Math.Exp(numerator / denominator);

            block.r_equivalent = r_equivalnt;
        }

        public static void setWI(BaseBlock block)
        {
            double K = 0;
            // starting from index = 2, to avoid top and bottom directions permeabilities.
            for (int i = 2; i < block.permeability_list.Length; i++)
            {
                K *= block.permeability_list[i] > 0 ? block.permeability_list[i] : 1;
            }
            K = Math.Sqrt(K);

            block.WI = 2 * Global.PI * K * block.h / (Math.Log(block.r_equivalent / block.well_radius) + block.skin);
        }

        public static double calculatePwf(BaseBlock block, double pressure, double Kr, double viscosity)
        {
            return pressure - (block.specified_flow_rate / (Kr / viscosity) * block.WI);
        }

        public static double calculateFlow_Rate(double P, double Pwf, double Kr, double viscosity, double WI)
        {
            return (P - Pwf) * (Kr / viscosity) * WI;
        }
    }
}
