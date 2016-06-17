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
                R[counter + 2] = -block.calculateR(data.grid, Global.Phase.Water, data.time_step, data.solubleGasPresent);

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
                jacobians[counter + 2] = new double[size];

                #region Oil
                // with respect to P
                jacobians[counter][data.phases.Length * block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter][data.phases.Length * block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                // with respect to Sw
                jacobians[counter][data.phases.Length * block.index + 2] = (block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, -1, Global.Variable.Saturation_Water, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                        // with respect to Sw
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.calculateDerivative(data.grid, Global.Phase.Oil, data.time_step, j, Global.Variable.Saturation_Water, data.solubleGasPresent) + minus_R[counter]) / Global.epsilon;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][data.phases.Length * block.index] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter + 1][data.phases.Length * block.index + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                // with respect to Sw
                jacobians[counter + 1][data.phases.Length * block.index + 2] = (block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, -1, Global.Variable.Saturation_Water, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j]] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = ( block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Saturation_Gas, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sw
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.calculateDerivative(data.grid, Global.Phase.Gas, data.time_step, j, Global.Variable.Saturation_Water, data.solubleGasPresent) + minus_R[counter + 1]) / Global.epsilon;
                    }
                }
                #endregion
                #region Water
                jacobians[counter + 2][data.phases.Length * block.index] = (block.calculateDerivative(data.grid, Global.Phase.Water, data.time_step, -1, Global.Variable.Pressure, data.solubleGasPresent) + minus_R[counter + 2]) / Global.epsilon;
                jacobians[counter + 2][data.phases.Length * block.index + 1] = 0;
                jacobians[counter + 2][data.phases.Length * block.index + 2] = (block.calculateDerivative(data.grid, Global.Phase.Water, data.time_step, -1, Global.Variable.Saturation_Water, data.solubleGasPresent) + minus_R[counter + 2]) / Global.epsilon;
                #endregion

                counter += data.phases.Length;
            }

            return jacobians;

        }

        public static double[] solveForDelta(double[][] jacobian, double[] minus_R)
        {
            int size = minus_R.Length;

            var vector_B = new DenseVector(minus_R);
            //var matrix_A = DenseMatrix.OfRowArrays(jacobian);
            var matrix_A = SparseMatrix.OfRowArrays(jacobian);

            var iteration_count_stop_criterion = new IterationCountStopCriterion<double>(1000);
            var residual_stop_criterion = new ResidualStopCriterion<double>(1e-10);
            var monitor = new Iterator<double>(iteration_count_stop_criterion, residual_stop_criterion);
            var solver = new MathNet.Numerics.LinearAlgebra.Double.Solvers.TFQMR();

            var delta = matrix_A.SolveIterative(vector_B, solver, monitor);

            //var delta = matrix_A.Solve(vector_B);

            return delta.ToArray();
        }

        public static void iterativeSolver(SimulationData data)
        {
            double[] minus_R;
            double[][] jacobian;
            double[] delta;

            double[] convergenceError = new double[2];

            resetTimeStep(data);

            //minus_R = calculate_minus_R(data);
            //jacobian = calculateJacobians(data, minus_R);

            do
            {

                minus_R = calculate_minus_R(data);
                jacobian = calculateJacobians(data, minus_R);
                delta = solveForDelta(jacobian, minus_R);

                convergenceError[0] = convergenceError[1];
                convergenceError[1] = checkTolerance(minus_R);

                bool repeat = stabilize_newton(delta, ref data.relaxation_factor, convergenceError, data);

                if (!repeat)
                {
                    updatePropertiesFromDelta(1, delta, data);
                }
                else
                {
                    for (int i = 0; i < data.grid.Length; i++)
                    {
                        data.grid[i].reset_n1(data.pvt, data.kr, data.porosity);
                    }

                    data.relaxation_factor = data.original_relaxation_factor;
                    convergenceError[1] = data.tolerance + 1;
                }


            } while (convergenceError[1] > data.tolerance);


            double MBE_Oil = MBE.checkOil(data);
            data.MBE_Gas = MBE.checkGas(data);
            double MBE_Water = MBE.checkWater(data);

            updateProperties(data);
        }

        private static void resetTimeStep(SimulationData data)
        {
            if (data.time_step > 0)
            {
                data.time_step = data.time_step * 2;
            }
            else
            {
                data.time_step = data.original_time_step;
            }
        }

        private static void updateProperties(SimulationData data)
        {
            double P, So, Sg, Sw;

            data.relaxation_factor = data.original_relaxation_factor;

            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1];
                Sg = data.grid[i].Sg[1];
                Sw = data.grid[i].Sw[1];
                So = 1 - Sw - Sg;

                data.grid[i].updateProperties(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);
            }
        }

        public static void updatePropertiesFromDelta(int time_level, double[] delta, SimulationData data)
        {
            double P, So, Sg, Sw;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1] + delta[counter];
                //Sg = data.grid[i].Sg[1] + delta[counter + 1] >= 0 ? data.grid[i].Sg[1] + delta[counter + 1] : 0;
                Sg = 0;
                Sw = data.grid[i].Sw[1] + delta[counter + 2] >= 0 ? data.grid[i].Sw[1] + delta[counter + 2] : 0;
                So = 1 - Sw - Sg;

                if (time_level == 0)
                {
                    data.grid[i].updateProperties(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);
                }
                else if (time_level == 1)
                {
                    data.grid[i].updateProperties_n1_k1(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);
                }

                counter += data.phases.Length;
            }
        }

        private static double checkTolerance(double[] delta)
        {
            double max = delta.Select(x => Math.Abs(x)).Max();
            return max;
        }

        private static bool stabilize_newton(double[] delta_x,ref double relaxation_factor, double[] convergenceError, SimulationData data)
        {
            if (convergenceError[0] != 0 && ((convergenceError[1] > convergenceError[0]) || convergenceError[1] / convergenceError[0] > 0.5))
            {
                relaxation_factor -= 0.05;

                if (relaxation_factor < data.minimum_relaxation)
                {
                    data.time_step *= 0.5;

                    return true;
                }

            }

            for (int i = 0; i < delta_x.Length; i++)
            {
                delta_x[i] = relaxation_factor * delta_x[i];
            }

            return false;
        }

        public static void RunSimulation(SimulationData data)
        {
            double end_time = 2 * 365;

            for (double current_time = 0; current_time < end_time; current_time += data.time_step)
            {
                iterativeSolver(data);

                Console.WriteLine(data.grid[299].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.grid[0].P[0] + ", " + data.MBE_Gas + ", " + data.time_step + ", " + data.grid[299].BHP[1]);
            }
        }

    }
}
