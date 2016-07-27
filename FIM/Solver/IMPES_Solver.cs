using FIM.Core;
using FIM.MaterialBalance;
using System;
using FIM.Extensions;
using FIM.Extensions.FullyImplicit;
using FIM.Extensions.IMPES;

using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FIM.Well;

namespace FIM.Solver
{
    class IMPES_Solver
    {
        // the entry point

        /// <summary>
        /// The entry point of the fully implicit simulation cycle.
        /// </summary>
        /// <param name="data">The <see cref="SimulationData"/>.</param>
        public static void RunSimulation(SimulationData data)
        {
            // size of the model grid.
            int size = data.grid.Length;

            // initialize the pressure coefficients matrix.
            double[][] pressureCoefficients = new double[size][];
            for (int i = 0; i < pressureCoefficients.Length; i++)
            {
                pressureCoefficients[i] = new double[size];
            }

            // initalize the empty constants matrix.
            double[] constants = new double[size];

            // intialize the empty delta p matrix. this matrix will store the solution of the Ax=B equation "where A is the pressure coefficients and B is the constant".
            double[] newP = new double[size];

            // the total simulated time.
            double end_time = data.endingTime;

            // the simulation loop.

            double current_time = 0;

            for (; current_time <= end_time;)
            {
                IterativeSolver(data, pressureCoefficients, constants, newP);
                // this way of updating the loop iterator "current_time" is used to correctly display current time
                // after the iterations. This is because the time step length may change during the iteration.
                current_time += data.timeStep;

                //Console.WriteLine(current_time + ", " + data.grid[0].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.grid[0].Rso[0]);
                Console.WriteLine(current_time + ", " + data.grid[0].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.MBE_Gas + ", " + data.wells[0].BHP[1] + ", " + data.wells[0].q_free_gas[0] + ", " + data.wells[0].q_water[0] + ", " + data.wells[0].q_oil[0]);
                //Console.WriteLine(data.MBE_Gas);
                //Console.WriteLine(data.grid[0].P[0] + ", " + data.wells[0].BHP[0]);
                Console.ReadKey();
            }
        }

        // the iteration cycle

        private static void IterativeSolver(SimulationData data, double[][] pressureCoefficients, double[] constants, double[] deltaP)
        {
            // an array containing previous "index 0" and current "index 1" material balance errors.
            double[] convergenceError = new double[2];

            // resets time step to its original value at the beginning of each new time step.
            ResetTimeStep(data);

            // non-linear iteration counter.
            int counter = 0;

            do
            {
                // formulation.
                IMPES_Preparation.CalculateConstantsMatrix(data, constants);
                IMPES_Preparation.CalculatePressureCoefficientsMatrix(data, pressureCoefficients);
                IMPES_Preparation.AddWellTerms(data, pressureCoefficients, constants);

                // solving the equations.
                deltaP = SolveForNewP(pressureCoefficients, constants);

                counter = Stabilize(data, convergenceError, deltaP, counter);

            } while (convergenceError[1] >= data.MBE_Tolerance && counter <= data.maximumNonLinearIterations);

            data.MBE_Oil = MBE.CheckOil(data);
            data.MBE_Gas = MBE.CheckGas(data);
            data.MBE_Water = MBE.CheckWater(data);

            UpdateProperties(data);
        }

        // solves the Ax=B equation "where A is the Jacobi and B is minus R" for deltas.
        private static double[] SolveForNewP(double[][] pressureCoefficients, double[] constants)
        {
            // note that this method while calculating delta_x "from the Jacobi and minus R matrices", it modifies them as a side effect.
            // this should be taken with much care.
            Mathematics.SolveLinearEquation.direct(pressureCoefficients, constants);
            return constants;
        }

        // contains the algorithms of convergence checking, relaxation of deltas and time cut off.
        private static int Stabilize(SimulationData data, double[] convergenceError, double[] delta, int counter)
        {
            UpdatePropertiesFromNewP(data, delta);

            // check convergence errors.
            bool repeat = CheckConvergence(data, convergenceError);
            //bool repeat = false;

            if (!repeat && (data.maximumNonLinearIterations - counter) != 1)
            {
                // increment the counter.
                counter += 1;
            }
            else
            {
                data.timeStep *= data.timeStepSlashingFactor;

                for (int i = 0; i < data.grid.Length; i++)
                {
                    data.grid[i].Reset(data);
                }

                // reset the convergenceError array so that it violates the convergence criteria
                convergenceError[0] = 0;
                convergenceError[1] = data.MBE_Tolerance;
                //convergenceError[1] = 0;

                // reset the counter.
                // this way we begin repeating the same time step with all the number of non linear iterations.
                counter = 0;
            }

            return counter;
        }

        // miscellaneous internal helper methods

        // returns the absolute value of the gas material balance error of the whole model.
        private static bool CheckConvergence(SimulationData data, double[] convergenceError)
        {
            bool firstIteration = convergenceError[0] == 0;

            convergenceError[0] = convergenceError[1];

            double tempMBE_Oil = Math.Abs(MBE.CheckOil(data));
            double tempMBE_Gas = Math.Abs(MBE.CheckGas(data));
            double tempMBE_Water = Math.Abs(MBE.CheckWater(data));

            //Math.Max(tempMBE_Gas, Math.Max(tempMBE_Oil, tempMBE_Water));

            convergenceError[1] = Math.Abs(MBE.CheckGas(data));

            bool MBE_Increasing = (convergenceError[1] > convergenceError[0]);
            //bool slowConvergence = convergenceError[1] / convergenceError[0] > data.maximumConvergenceErrorRatio;

            // if this is not the first iteration and either one of the following conditions is true:
            // 1- material balance error is increasing not decreasing.
            // 2- the rate of convergence is slow.
            if (!firstIteration && (MBE_Increasing /*|| slowConvergence*/))
            {
                return true;
            }

            return false;
        }

        // resets the time step to the original value.
        // this code may be modified to make gradual return from the cut off time step to the original value.
        private static void ResetTimeStep(SimulationData data)
        {
            data.timeStep = data.originalTimeStep;
        }

        // updates the n1 time level properties with better approximation after solving the Ax=B equation.
        private static void UpdatePropertiesFromNewP(SimulationData data, double[] newP)
        {
            BaseBlock block, neighbor_block, upstream_block, downstream_block;
            double P, So, Sg, Sw;
            double kr, B, viscosity, RsO;
            double temp, transmissibility;
            double q_oil = 0, q_free_gas = 0, q_water = 0;

            int neighbor_block_index;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                P = newP[i];
                // assert p is always greater than 0.
                P = P > 0 ? P : 0;


                if (block.type == Global.BlockType.WellBlock)
                {
                    q_oil = data.wells[i].q_oil[0];
                    q_free_gas = data.wells[i].q_free_gas[1];
                    q_water = data.wells[i].q_water[0];
                }
                else
                {
                    q_oil = 0;
                    q_free_gas = 0;
                    q_water = 0;
                }

                // oil phase
                temp = 0;
                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    neighbor_block_index = block.neighbourBlocksIndices[j];
                    if (neighbor_block_index < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[neighbor_block_index];
                    transmissibility = block.transmissibility_list[j];

                    if (block.P[1] >= neighbor_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    kr = upstream_block.Kro[0];
                    B = 0.5 * (block.Bo[0] + neighbor_block.Bo[0]);
                    viscosity = 0.5 * (block.viscosityOil[0] + neighbor_block.viscosityOil[0]);

                    temp += transmissibility * kr / (viscosity * B) * (newP[neighbor_block_index] - newP[i]);
                }
                So = block.Bo[1] / block.Vp[1] * (block.Vp[0] * block.So[0] / block.Bo[0] + data.timeStep * Global.a * (temp - q_oil));
                // assert So is always greater than 0.
                So = So > 1 ? 1 : So;

                // water phase
                temp = 0;
                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    neighbor_block_index = block.neighbourBlocksIndices[j];
                    if (neighbor_block_index < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[neighbor_block_index];
                    transmissibility = block.transmissibility_list[j];

                    if (block.P[1] >= neighbor_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    kr = upstream_block.Krw[0];
                    B = 0.5 * (block.Bw[0] + neighbor_block.Bw[0]);
                    viscosity = 0.5 * (block.viscosityWater[0] + neighbor_block.viscosityWater[0]);

                    temp += transmissibility * kr / (viscosity * B) * (newP[neighbor_block_index] - newP[i]);
                }
                Sw = block.Bw[1] / block.Vp[1] * (block.Vp[0] * block.Sw[0] / block.Bw[0] + data.timeStep * Global.a * (temp - q_water));
                // assert Sw is always greater than 0.
                Sw = Sw > 1 ? 1 : Sw;

                temp = 0;
                // free gas
                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    neighbor_block_index = block.neighbourBlocksIndices[j];
                    if (neighbor_block_index < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[neighbor_block_index];
                    transmissibility = block.transmissibility_list[j];

                    if (block.P[1] >= neighbor_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    kr = upstream_block.Krg[0];
                    B = 0.5 * (block.Bg[0] + neighbor_block.Bg[0]);
                    viscosity = 0.5 * (block.viscosityGas[0] + neighbor_block.viscosityGas[0]);

                    temp += transmissibility * kr / (viscosity * B) * (newP[neighbor_block_index] - newP[i]);
                }

                // solution gas
                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    neighbor_block_index = block.neighbourBlocksIndices[j];
                    if (neighbor_block_index < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[neighbor_block_index];
                    transmissibility = block.transmissibility_list[j];

                    if (block.P[1] >= neighbor_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    kr = upstream_block.Kro[0];
                    B = 0.5 * (block.Bo[0] + neighbor_block.Bo[0]);
                    viscosity = 0.5 * (block.viscosityOil[0] + neighbor_block.viscosityOil[0]);
                    RsO = 0.5 * (block.Rso[0] + neighbor_block.Rso[0]);

                    temp += RsO * transmissibility * kr / (viscosity * B) * (newP[neighbor_block_index] - newP[i]);

                    //temp += -1 * block.Rso[1] * transmissibility * kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1]);
                }

                double term1 = block.Vp[0] * block.Sg[0] / block.Bg[0];
                double term2 = block.Vp[0] * block.So[0] / block.Bo[0] * block.Rso[0];
                double term3 = block.Vp[1] * block.So[1] / block.Bo[1] * block.Rso[1];
                double term4 = data.timeStep * Global.a * (temp - q_free_gas - q_oil * block.Rso[1]);
                Sg = block.Bg[1] / block.Vp[1] * (term1 + term2 - term3 + term4);
                
                Sg = Sg > 0 ? Sg : 0;

                double s = Sg + Sw + So;

                block.UpdateProperties(data, P, Sw, Sg, So, 1);
            }

        }

        // updates the new time step properties after convergence.
        private static void UpdateProperties(SimulationData data)
        {
            double P, So, Sg, Sw;

            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1];
                Sg = data.grid[i].Sg[1];
                Sw = data.grid[i].Sw[1];
                So = 1 - Sw - Sg;

                data.grid[i].UpdateProperties(data, P, Sw, Sg, So, 0);
            }
        }
    }
}
