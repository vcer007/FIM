using FIM.Core;
using FIM.MaterialBalance;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Solver
{
    class FIM_block
    {
        public static double[] solveForDelta(double[][] jacobian, double[] minus_R)
        {
            double[] delta = new double[minus_R.Length];
            for (int i = 0; i < minus_R.Length; i++)
            {
                delta[i] = minus_R[i];
            }

            Mathematics.SolveLinearEquation.direct(jacobian, delta);
            return delta;
        }

        public static void iterativeSolver(SimulationData data, double[][] jacobian, double[] minus_R, double[] delta)
        {
            double[] convergenceError = new double[2];

            resetTimeStep(data);


            int counter = 0;
            do
            {

                minus_R = Extensions.SingleBlockExtensions.calculate_minus_R(data);
                jacobian = Extensions.SingleBlockExtensions.calculateJacobians(data, minus_R);

                //Console.Write("Max minus_R : " + minus_R.Select(x => Math.Abs(x)).Max() + ", ");

                delta = solveForDelta(jacobian, minus_R);

                convergenceError[0] = convergenceError[1];

                convergenceError[1] = checkTolerance(data);

                bool repeat = stabilize_newton(delta, ref data.relaxation_factor, convergenceError, data);
                //bool repeat = false;

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
                Sw = data.grid[i].Sw[0];
                So = 1 - Sw - Sg;

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
                //Console.WriteLine("Sg : " + Sg + " P : " + P + " MBG : " + MBE.checkGas(data));

                Sg = Sg > 0 ? Sg : 0;
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

            double end_time = 1 * 178;

            for (double current_time = 0; current_time <= end_time; current_time += data.time_step)
            {
                iterativeSolver(data, jacobian, minus_R, delta);



                //Console.WriteLine("###################################################################");
                Console.WriteLine(current_time + ", " + data.grid[0].P[0] + ", " + data.grid[0].Rso[0] + ", " + data.grid[0].Sg[0] + ", " + data.MBE_Gas + ", " + data.grid[0].GOR + ", " + data.grid[0].BHP[0]);
                //Console.WriteLine("###################################################################");
                //Console.ReadKey();
            }
        }
    }
}
