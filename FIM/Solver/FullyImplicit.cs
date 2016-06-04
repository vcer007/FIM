using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Solvers;
using FIM.Core;
using FIM.Rock;
using FIM.Fluid;
using FIM.MaterialBalance;

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

                R[counter] = -block.calculateR(data.grid, Global.Phase.Oil, data.time_step, data.solubleGasPresent);
                R[counter + 1] = -block.calculateR(data.grid, Global.Phase.Gas, data.time_step, data.solubleGasPresent);

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
                jacobians[counter][data.phases.Length * block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                // with respect to Sg

                jacobians[counter][data.phases.Length * block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][data.phases.Length * block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter + 1][data.phases.Length * block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
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

            var vector_B = new DenseVector(minus_R);
            var matrix_A = DenseMatrix.OfRowArrays(jacobian);

            //var iteration_count_stop_criterion = new IterationCountStopCriterion<double>(1000);
            //var residual_stop_criterion = new ResidualStopCriterion<double>(1e-10);
            //var monitor = new Iterator<double>(iteration_count_stop_criterion, residual_stop_criterion);
            //var solver = new TFQMR();

            //var delta = matrix_A.SolveIterative(vector_B, solver, monitor);

            var delta = matrix_A.Solve(vector_B);

            return delta.ToArray(); ;
        }

        public static void iterativeSolver(SimulationData data)
        {
            double[] minus_R;
            double[][] jacobian;
            double[] delta;

            minus_R = calculate_minus_R(data);
            jacobian = calculateJacobians(data, minus_R);

            do
            {

                minus_R = calculate_minus_R(data);
                delta = solveForDelta(jacobian, minus_R);
                updatePropertiesFromDelta(1, delta, data);

            } while (checkTolerance(minus_R) > data.tolerance);


            double MBE_Oil = MBE.checkOil(data);
            double MBE_Gas = MBE.checkGas(data);

            updateProperties(data);
        }

        private static void updateProperties(SimulationData data)
        {
            double P, So, Sg, Sw;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1];
                Sg = data.grid[i].Sg[1];
                So = 1 - data.grid[i].Sw[1] - Sg;
                Sw = data.grid[i].Sw[1];

                data.grid[i].updateProperties(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);

                counter += 2;
            }
        }

        public static void updatePropertiesFromDelta(int time_level, double[] delta, SimulationData data)
        {
            double P, So, Sg, Sw;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1] + 1 * delta[counter];
                Sg = data.grid[i].Sg[1] + 1 * delta[counter + 1] >= 0 ? data.grid[i].Sg[1] + 1 * delta[counter + 1] : 0;
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

                counter += 2;
            }
        }

        private static double checkTolerance(double[] delta)
        {
            return delta.Select(x => Math.Abs(x)).Max();
        }

        public static void RunSimulation(SimulationData data)
        {
            double end_time = 2 * 365;

            for (double current_time = 0; current_time < end_time; current_time += data.time_step)
            {
                iterativeSolver(data);

                Console.WriteLine(data.grid[299].P[0] + ", " + data.grid[299].Sg[0] + ", " + data.grid[299].BHP[0]);
            }
        }

    }
}
