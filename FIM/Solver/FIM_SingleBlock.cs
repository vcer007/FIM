using FIM.Core;
using FIM.Extensions;
using FIM.MaterialBalance;
using MathNet.Numerics.LinearAlgebra.Complex;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Solver
{
    class FIM_SingleBlock
    {
        public static void calculate_minus_R(SimulationData data, double[] minus_R)
        {
            //int size = data.x * data.y * data.z * data.phases.Length;
            //double[] minus_R = new double[size];

            BaseBlock block;
            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                minus_R[counter] = -1 * block.single_block_R(data, Global.Phase.Oil);
                minus_R[counter + 1] = -1 * block.single_block_R(data, Global.Phase.Gas);
                //minus_R[counter + 2] = -1 * block.single_block_R(data, Global.Phase.Water);

                counter += data.phases.Length;
            }
        }

        public static double[][] calculateJacobians(SimulationData data, double[] minus_R, double[][] jacobian)
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
                //jacobians[counter + 2] = new double[size];

                #region Oil
                // with respect to P
                jacobians[counter][data.phases.Length * block.index] = block.qasim_jacobian(data, Global.Phase.Oil, -1, Global.Variable.Pressure);
                // with respect to Sg
                jacobians[counter][data.phases.Length * block.index + 1] = block.qasim_jacobian(data, Global.Phase.Oil, -1, Global.Variable.Saturation_Gas);
                // with respect to Sw
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][data.phases.Length * block.index] = block.qasim_jacobian(data, Global.Phase.Gas, -1, Global.Variable.Pressure);
                // with respect to Sg
                jacobians[counter + 1][data.phases.Length * block.index + 1] = block.qasim_jacobian(data, Global.Phase.Gas, -1, Global.Variable.Saturation_Gas);

                #endregion

                counter += data.phases.Length;
            }

            return jacobians;

        }

        public static double[] solveForDelta(double[][] jacobian, double[] minus_R)
        {
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

                calculate_minus_R(data, minus_R);
                jacobian = calculateJacobians(data, minus_R, jacobian);

                Console.Write("Max minus_R : " + minus_R.Select(x => Math.Abs(x)).Max() + ", ");

                delta = solveForDelta(jacobian, minus_R);

                convergenceError[0] = convergenceError[1];

                convergenceError[1] = checkTolerance(data);

                //bool repeat = stabilize_newton(delta, ref data.relaxation_factor, convergenceError, data);
                bool repeat = false;

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

                counter += 1;
            } while (convergenceError[1] > data.tolerance && counter <= 25);
            // 

            data.MBE_Oil = MBE.checkOil(data);
            data.MBE_Gas = MBE.checkGas(data);

            updateProperties(data);
        }

        private static void resetTimeStep(SimulationData data)
        {
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
                data.grid[i].updateProperties(data.pvt, data.kr, data.porosity, P, Sw, So, Sg, p_increment);
            }
        }

        public static void updatePropertiesFromDelta(int time_level, double[] delta, SimulationData data)
        {
            double P, So, Sg, Sw;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                //P = data.grid[i].P[1] + delta[counter] > 0 ? data.grid[i].P[1] + delta[counter] : 0;
                P = data.grid[i].P[1] + delta[counter];

                //if (P <= data.pvt.bubble_point_pressure)
                //{
                //    data.grid[i].saturated = true;
                //}
                //if (data.grid[i].saturated)
                //{
                //    if (P > data.pvt.bubble_point_pressure)
                //    {
                //        P = data.pvt.bubble_point_pressure;
                //    }
                //}

                Sg = data.grid[i].Sg[1] + delta[counter + 1];
                Console.WriteLine("Sg : " + Sg + " P : " + P);

                Sg = data.grid[i].Sg[1] + delta[counter + 1] > 0 ? data.grid[i].Sg[1] + delta[counter + 1] : data.grid[i].Sg[1];
                //Sw = data.grid[i].Sw[1] + delta[counter + 2];
                Sw = 0;

                So = 1 - Sw - Sg;

                data.grid[i].updateProperties_n1_k1(data.pvt, data.kr, data.porosity, P, Sw, So, Sg);

                counter += data.phases.Length;
            }

        }

        private static double checkTolerance(SimulationData data)
        {
            double temmp = Math.Abs(MBE.checkGas(data));
            //Console.Write("tolerance MBG : " + temmp + ", " + "P : " + data.grid[0].P[1]);
            //Console.WriteLine();
            return temmp;
        }

        private static bool stabilize_newton(double[] delta_x, ref double relaxation_factor, double[] convergenceError, SimulationData data)
        {
            if (convergenceError[0] != 0 && ((convergenceError[1] > convergenceError[0]) || convergenceError[1] / convergenceError[0] > 0.5))
            {
                relaxation_factor -= 0.1;
                Console.WriteLine("relax : " + relaxation_factor);

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

            double end_time = 8 * 365;

            for (double current_time = 0; current_time < end_time; current_time += data.time_step)
            {
                iterativeSolver(data, jacobian, minus_R, delta);



                Console.WriteLine("###################################################################");
                Console.WriteLine(current_time + ", " + data.grid[0].P[0] + ", " + data.grid[0].Rso[0] + ", " + data.grid[0].Sg[0] + ", " + data.MBE_Gas + ", " + data.time_step + ", " + data.grid[0].BHP[1]);
                Console.WriteLine("###################################################################");
                Console.ReadKey();
            }
        }
    }
}
