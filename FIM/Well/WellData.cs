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

            double distance = 0;
            double numerator = 0, denominator = 0;

            BaseBlock neighbour;
            // neighbour block internal index of "block".
            int index;

            for (int i = 2; i < block.neighbour_blocks_indices.Length; i++)
            {
                if (block.neighbour_blocks_indices[i] >= 0)
                {
                    neighbour = grid[block.neighbour_blocks_indices[i]];
                    index = Array.IndexOf(neighbour.neighbour_blocks_indices, block.index);

                    distance = (block.delta_x_list[i] + neighbour.delta_x_list[index]) / 2;
                    numerator += block.boundary_length_list[i] / distance * Math.Log(distance) - 0.5 * Math.PI;
                    denominator += (block.boundary_length_list[i]) / distance;
                }
            }

            r_equivalnt = Math.Exp(numerator / denominator);

            if (grid.Length == 1)
            {
                r_equivalnt = 0.2 * block.delta_x_list[3];
            }

            block.r_equivalent = r_equivalnt;
        }

        public static void setWI(BaseBlock block)
        {
            double K = 0;
            // starting from index = 2, to avoid top and bottom directions permeabilities.
            int counter = 0;
            for (int i = 2; i < block.permeability_list.Length; i++)
            {
                if (block.permeability_list[i] > 0)
                {
                    K += block.permeability_list[i];
                    counter += 1;
                }
                else
                {
                    
                }
            }

            K = K / counter;

            block.WI = 2 * Global.Bc * Global.PI * K * block.h / (Math.Log(block.r_equivalent / block.well_radius) + block.skin);
        }

        public static double calculatePwf(BaseBlock block, double pressure, double Kr, double viscosity, double FVF)
        {
            double mobility = Kr / (viscosity * FVF);
            return pressure - ( block.q_oil[0] / (mobility * block.WI) );
        }

        public static double calculateFlow_Rate(double P, double Pwf, double Kr, double viscosity, double WI, double FVF)
        {
            return (P - Pwf) * (Kr / (viscosity * FVF)) * WI;
            //return 0;
        }
    }
}
