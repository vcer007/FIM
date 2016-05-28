using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FIM.Core;

namespace FIM.Solver
{
    class FullyImplicit
    {
        public static double[][] calculateJacobians(SimulationData data, double[] R)
        {
            int size = data.x * data.y * data.z * data.phases.Length;
            double[][] jacobians = new double[size][];

            BaseBlock block;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];
                jacobians[counter] = new double[size];
                jacobians[counter + 1] = new double[size];
                jacobians[counter + 2] = new double[size];

                // oil
                jacobians[counter][block.index] = block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_stpe, -1, Global.Variable.Pressure, true);
                jacobians[counter][block.index + 1] = block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_stpe, -1, Global.Variable.Saturation, true);
                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        jacobians[counter][3 * block.neighbour_blocks_indices[j]] = block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_stpe, j, Global.Variable.Pressure, true);
                        jacobians[counter][3 * block.neighbour_blocks_indices[j] + 1] = block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_stpe, j, Global.Variable.Saturation, true);
                    }
                }
                // gas
                jacobians[counter + 1][block.index] = block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_stpe, -1, Global.Variable.Pressure, true);
                jacobians[counter + 1][block.index + 2] = block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_stpe, -1, Global.Variable.Pressure, true);
                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        jacobians[counter + 1][3 * block.neighbour_blocks_indices[j]] = block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_stpe, j, Global.Variable.Pressure, true);
                        jacobians[counter + 1][3 * block.neighbour_blocks_indices[j] + 1] = block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_stpe, j, Global.Variable.Saturation, true);
                    }
                }
                // water
                jacobians[counter + 2][block.index] = block.calculateDerivative(data.grid, Global.Phase.Water, data.time_stpe, -1, Global.Variable.Pressure, true);
                jacobians[counter + 2][block.index + 1] = block.calculateDerivative(data.grid, Global.Phase.Water, data.time_stpe, -1, Global.Variable.Saturation, true);
                jacobians[counter + 2][block.index + 2] = block.calculateDerivative(data.grid, Global.Phase.Water, data.time_stpe, -1, Global.Variable.Saturation, true);

                counter += 3;
            }

            return jacobians;

        }

        public static double[] calculate_R(SimulationData data)
        {
            int size = data.x * data.y * data.z * data.phases.Length;
            double[] R = new double[size];

            BaseBlock block;
            int counter = 0;
            for (int i = 0; i < data.grid.Length; i ++)
            {
                block = data.grid[i];

                R[counter] = block.calculateR(data.grid, Global.Phase.Oil, data.time_stpe, true);
                R[counter + 1] = block.calculateR(data.grid, Global.Phase.Water, data.time_stpe, true);
                R[counter + 2] = block.calculateR(data.grid, Global.Phase.Gas, data.time_stpe, true);

                counter += 3;
            }

            return R;
        }
    }
}
