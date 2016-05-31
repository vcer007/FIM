﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FIM.Core;
using FIM.Rock;
using FIM.Fluid;

namespace FIM.Solver
{
    class FullyImplicit
    {

        public static double[] calculate_minus_R(SimulationData data)
        {
            int size = data.x * data.y * data.z * data.phases.Length;
            double[] R = new double[size];

            BaseBlock block;
            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                R[counter] = -block.calculateR(data.grid, Global.Phase.Oil, data.time_step, true);
                R[counter + 1] = -block.calculateR(data.grid, Global.Phase.Water, data.time_step, true);

                counter += data.phases.Length;
            }

            return R;
        }

        public static double[][] calculateJacobians(SimulationData data, double[] minus_R)
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
                jacobians[counter][block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Pressure, true) + minus_R[counter]) / Global.epsilon;
                // with respect to Sg

                jacobians[counter][block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Saturation_Gas, true) + minus_R[counter]) / Global.epsilon;
                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Pressure, true) + minus_R[counter]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Saturation_Gas, true) + minus_R[counter]) / Global.epsilon;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Pressure, true) + minus_R[counter + 1]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter + 1][block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Saturation_Gas, true) + minus_R[counter + 1]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Pressure, true) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Saturation_Gas, true) + minus_R[counter + 1]) / Global.epsilon;
                    }
                }
                #endregion

                counter += data.phases.Length;
            }

            return jacobians;

        }

        public static double[] solveForDelta(double[][] jacobian, double[] minus_R)
        {
            int size = minus_R.Length;
            double[] delta = new double[size];

            return delta;
        }

        public static void iterativeSolver(SimulationData data)
        {
            double[] minus_R;
            double[][] jacobian;
            double[] delta;

            minus_R = calculate_minus_R(data);
            jacobian = calculateJacobians(data, minus_R);
            delta = solveForDelta(jacobian, minus_R);
            updatePropertiesFromDelta(1, delta, data);

            while (checkTolerance(delta) > data.tolerance)
            {
                minus_R = calculate_minus_R(data);
                jacobian = calculateJacobians(data, minus_R);
                delta = solveForDelta(jacobian, minus_R);
                updatePropertiesFromDelta(1, delta, data);
            }

            updatePropertiesFromDelta(0, delta, data);

        }

        public static void updatePropertiesFromDelta(int time_level, double[] delta, SimulationData data)
        {
            double P, So, Sg, Sw;

            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1] + delta[1];
                Sg = data.grid[i].Sg[1] + delta[2];
                So = 1 - data.grid[i].Sw[1] - Sg;
                Sw = data.grid[i].Sw[1];

                if (time_level == 0)
                {
                    data.grid[i].updateProperties(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);
                }
                else if (time_level == 1)
                {
                    data.grid[i].updateProperties_n1(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);
                }
            }
        }

        private static double checkTolerance(double[] delta)
        {
            return delta.Min();
        }

    }
}
