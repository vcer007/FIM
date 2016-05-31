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

                #region Oil
                // with respect to P
                jacobians[counter][block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Pressure, true) - R[counter]) / Global.epsilon;
                // with respect to Sg

                jacobians[counter][block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Saturation_Gas, true) - R[counter]) / Global.epsilon;
                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Pressure, true) - R[counter]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Saturation_Gas, true) - R[counter]) / Global.epsilon;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Pressure, true) - R[counter + 1]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter + 1][block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Saturation_Gas, true) - R[counter + 1]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Pressure, true) - R[counter + 1]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Saturation_Gas, true) - R[counter + 1]) / Global.epsilon;
                    }
                }
                #endregion

                counter += data.phases.Length;
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

                R[counter] = block.calculateR(data.grid, Global.Phase.Oil, data.time_step, true);
                R[counter + 1] = block.calculateR(data.grid, Global.Phase.Water, data.time_step, true);

                counter += data.phases.Length;
            }

            return R;
        }
    }
}
