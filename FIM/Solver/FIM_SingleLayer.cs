using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Solvers;
using FIM.Core;
using FIM.RockData;
using FIM.FluidData;
using FIM.MaterialBalance;
using FIM.Extensions;

namespace FIM.Solver
{
    class FIM_SingleLayer
    {

        //public static void calculate_minus_R(SimulationData data, double[] minus_R)
        //{
        //    //int size = data.x * data.y * data.z * data.phases.Length;
        //    //double[] minus_R = new double[size];

        //    BaseBlock block;
        //    int counter = 0;
        //    for (int i = 0; i < data.grid.Length; i++)
        //    {
        //        block = data.grid[i];

        //        double temp = -1 * block.calculateR(data.grid, Global.Phase.Oil, data.time_step, data.solubleGasPresent);
        //        minus_R[counter] = temp;
        //        minus_R[counter + 1] = -1 * block.calculateR(data.grid, Global.Phase.Gas, data.time_step, data.solubleGasPresent);
        //        minus_R[counter + 2] = -1 * block.calculateR(data.grid, Global.Phase.Water, data.time_step, data.solubleGasPresent);

        //        counter += data.phases.Length;
        //    }
        //}

        //public static double[][] calculateJacobians(SimulationData data, double[] minus_R, double[][] jacobian)
        //{
        //    int size = data.grid.Length * data.phases.Length;
        //    double[][] jacobians = new double[size][];

        //    BaseBlock block;

        //    int counter = 0;
        //    for (int i = 0; i < data.grid.Length; i++)
        //    {
        //        block = data.grid[i];

        //        jacobians[counter] = new double[size];
        //        jacobians[counter + 1] = new double[size];
        //        jacobians[counter + 2] = new double[size];

        //        #region Oil
        //        // with respect to P
        //        jacobians[counter][data.phases.Length * block.index] = (block.perturb(data, Global.Phase.Oil, -1, Global.Variable.Pressure) + minus_R[counter]) / block.dp;
        //        // with respect to Sg
        //        jacobians[counter][data.phases.Length * block.index + 1] = (block.perturb(data, Global.Phase.Oil, -1, Global.Variable.Saturation_Gas) + minus_R[counter]) / block.dsg;
        //        // with respect to Sw
        //        jacobians[counter][data.phases.Length * block.index + 2] = (block.perturb(data, Global.Phase.Oil, -1, Global.Variable.Saturation_Water) + minus_R[counter]) / block.dsw;

        //        for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
        //        {
        //            if (block.neighbour_blocks_indices[j] >= 0)
        //            {
        //                // with respect to P
        //                jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j]] = (block.perturb(data, Global.Phase.Oil, j, Global.Variable.Pressure) + minus_R[counter]) / block.dp;
        //                // with respect to Sg
        //                jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = (block.perturb(data, Global.Phase.Oil, j, Global.Variable.Saturation_Gas) + minus_R[counter]) / block.dsg;
        //                // with respect to Sw
        //                jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.perturb(data, Global.Phase.Oil, j, Global.Variable.Saturation_Water) + minus_R[counter]) / block.dsw;
        //            }
        //        }
        //        #endregion
        //        #region Gas
        //        // with respect to P
        //        jacobians[counter + 1][data.phases.Length * block.index] = (block.perturb(data, Global.Phase.Gas, -1, Global.Variable.Pressure) + minus_R[counter + 1]) / block.dp;
        //        // with respect to Sg
        //        jacobians[counter + 1][data.phases.Length * block.index + 1] = (block.perturb(data, Global.Phase.Gas, -1, Global.Variable.Saturation_Gas) + minus_R[counter + 1]) / block.dsg;
        //        // with respect to Sw
        //        jacobians[counter + 1][data.phases.Length * block.index + 2] = (block.perturb(data, Global.Phase.Gas, -1, Global.Variable.Saturation_Water) + minus_R[counter + 1]) / block.dsw;

        //        for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
        //        {
        //            if (block.neighbour_blocks_indices[j] >= 0)
        //            {
        //                // with respect to P
        //                jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j]] = (block.perturb(data, Global.Phase.Gas, j, Global.Variable.Pressure) + minus_R[counter + 1]) / block.dp;
        //                // with respect to Sg
        //                jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = (block.perturb(data, Global.Phase.Gas, j, Global.Variable.Saturation_Gas) + minus_R[counter + 1]) / block.dsg;
        //                // with respect to Sw
        //                jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.perturb(data, Global.Phase.Gas, j, Global.Variable.Saturation_Water) + minus_R[counter + 1]) / block.dsw;
        //            }
        //        }
        //        #endregion
        //        #region Water
        //        jacobians[counter + 2][data.phases.Length * block.index] = (block.perturb(data, Global.Phase.Water, -1, Global.Variable.Pressure) + minus_R[counter + 2]) / block.dp;
        //        jacobians[counter + 2][data.phases.Length * block.index + 1] = 0;
        //        jacobians[counter + 2][data.phases.Length * block.index + 2] = (block.perturb(data, Global.Phase.Water, -1, Global.Variable.Saturation_Water) + minus_R[counter + 2]) / block.dsw;
        //        #endregion

        //        counter += data.phases.Length;
        //    }

        //    return jacobians;

        //}


        public static double[] solveForDelta(double[][] jacobian, double[] minus_R)
        {
            //int size = minus_R.Length;

            var vector_B = new DenseVector(minus_R);
            //var matrix_A = DenseMatrix.OfRowArrays(jacobian);
            var matrix_A = SparseMatrix.OfRowArrays(jacobian);

            //var iteration_count_stop_criterion = new IterationCountStopCriterion<double>(1000);
            //var residual_stop_criterion = new ResidualStopCriterion<double>(1e-10);
            //var monitor = new Iterator<double>(iteration_count_stop_criterion, residual_stop_criterion);
            //var solver = new MathNet.Numerics.LinearAlgebra.Double.Solvers.TFQMR();

            //var delta = matrix_A.SolveIterative(vector_B, solver, monitor);

            //var delta = matrix_A.Solve(vector_B);
            //return delta.ToArray();

            Mathematics.SolveLinearEquation.direct(jacobian, minus_R);
            return minus_R;
        }

        public static void iterativeSolver(SimulationData data, double[][] jacobian, double[] minus_R, double[] delta)
        {
            double[] convergenceError = new double[2];

            resetTimeStep(data);


            int counter = 0;
            do
            {
                //Console.WriteLine("GIP = " + MBE.GIP(data, 0));

                //calculate_minus_R(data, minus_R);
                //jacobian = calculateJacobians(data, minus_R, jacobian);

                //Console.WriteLine("Max minus_R : " + minus_R.Max());
                //Console.WriteLine("FGIP : " + MBE.FreeGasIP(data) + "   SGIP : " + MBE.SolubleGasIP(data) + "    Rso : " + data.grid[299].Rso[1]);

                NumericalPerturbation.calculateMinusR_Matrix(data, minus_R);
                NumericalPerturbation.calculateJacobi_Matrix(data, minus_R, jacobian);

                WellTerms.add(data, jacobian, minus_R);

                delta = solveForDelta(jacobian, minus_R);


                convergenceError[0] = convergenceError[1];

                convergenceError[1] = checkTolerance(data);

                bool repeat = stabilize_newton(delta, ref data.relaxation_factor, convergenceError, data);

                if (!repeat)
                {
                    updatePropertiesFromDelta(1, delta, data);
                }
                else
                {
                    for (int i = 0; i < data.grid.Length; i++)
                    {
                        data.grid[i].reset_n1(data);
                    }

                    data.relaxation_factor = data.original_relaxation_factor;
                    convergenceError[1] = data.tolerance + 1;
                }

                counter += 1;
            } while (convergenceError[1] > data.tolerance && counter <= 25);
            // 

            data.MBE_Oil = MBE.checkOil(data);
            data.MBE_Gas = MBE.checkGas(data);
            double MBE_Water = MBE.checkWater(data);

            updateProperties(data);
        }

        private static void resetTimeStep(SimulationData data)
        {
            //if (data.time_step > 0)
            //{
            //    data.time_step = data.time_step * 2;
            //}
            //else
            //{
            //    data.time_step = data.original_time_step;
            //}
            data.time_step = data.original_time_step;
        }

        private static void updateProperties(SimulationData data)
        {
            double P, So, Sg, Sw;
            double p_increment = 0;

            data.relaxation_factor = data.original_relaxation_factor;

            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1];
                Sg = data.grid[i].Sg[1];
                Sw = data.grid[i].Sw[1];
                So = 1 - Sw - Sg;

                //p_increment = data.grid[i].P[1] - data.grid[i].P[0];
                //p_increment = data.grid.Select(x => x.P[1] - x.P[0]).Max();
                //p_increment = -400;
                data.grid[i].updateProperties(data, P, Sw, Sg, 0);
            }
        }

        public static void updatePropertiesFromDelta(int time_level, double[] delta, SimulationData data)
        {
            double P, So, Sg, Sw;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                //P = data.grid[i].P[1] + delta[counter] < data.grid[i].P[0]  ? data.grid[i].P[1] + delta[counter] : data.grid[i].P[1];
                //P = P > 0 ? P : 0;
                P = data.grid[i].P[1] + delta[counter] > 0 ? data.grid[i].P[1] + delta[counter] : 0;

                Sg = data.grid[i].Sg[1] + delta[counter + 1] > 0 ? data.grid[i].Sg[1] + delta[counter + 1] : 0;

                //Sg = 0;
                Sw = data.grid[i].Sw[1] + delta[counter + 2];

                So = 1 - Sw - Sg;

                if (time_level == 0)
                {
                    data.grid[i].updateProperties(data, P, Sw, Sg, 0);
                }
                else if (time_level == 1)
                {
                    data.grid[i].updateProperties(data, P, Sw, Sg, 1);
                }

                counter += data.phases.Length;
            }

        }

        private static double checkTolerance(SimulationData data)
        {
            double temmp = Math.Abs(MBE.checkGas(data));
            //Console.WriteLine("tolerance MBG : " + temmp);
            return temmp;
        }

        private static bool stabilize_newton(double[] delta_x, ref double relaxation_factor, double[] convergenceError, SimulationData data)
        {
            if (convergenceError[0] != 0 && ((convergenceError[1] > convergenceError[0]) || convergenceError[1] / convergenceError[0] > 0.5))
            {
                relaxation_factor -= 0.1;
                //Console.WriteLine("relax : " + relaxation_factor);

                if (relaxation_factor < data.minimum_relaxation)
                {
                    data.time_step *= 0.5;
                    return true;
                }

            }

            if (relaxation_factor != 1)
            {
                for (int i = 0; i < delta_x.Length; i++)
                {
                    delta_x[i] = relaxation_factor * delta_x[i];
                }
            }

            return false;
        }

        public static void RunSimulation(SimulationData data)
        {
            int size = data.grid.Length * data.phases.Length;
            double[][] jacobian = new double[size][];
            for (int i = 0; i < jacobian.Length; i++)
            {
                jacobian[i] = new double[size];
            }
            double[] minus_R = new double[jacobian.Length];
            double[] delta = new double[jacobian.Length];

            double end_time = 5 * 46;

            for (double current_time = data.time_step; current_time <= end_time; current_time += data.time_step)
            {
                double temp = MBE.GIP(data, 0);

                iterativeSolver(data, jacobian, minus_R, delta);


                if (current_time == 360)
                {
                    //data.time_step = 60;
                }

                //Console.WriteLine("###################################################################");
                Console.WriteLine(current_time + ", " + data.grid[0].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.grid[0].Rso[0]);
                //Console.WriteLine(current_time + ", " + data.grid[0].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.grid[0].Rso[0] + ", " + data.wells[0].BHP[1] + ", " + data.wells[0].q_free_gas[0] + ", " + data.wells[0].q_solution_gas[0] + ", " + data.wells[0].q_oil[0]);
                //Console.WriteLine(data.grid[0].P[0] + ", " + data.wells[0].BHP[0]);
                Console.ReadKey();
                //Console.WriteLine("###################################################################");
            }
        }



    }
}
