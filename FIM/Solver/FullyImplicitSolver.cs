using FIM.Core;
using FIM.MaterialBalance;
using FIM.Extensions.FullyImplicit;
using FIM.Extensions;

using System;

/// <summary>
/// This name space contains the fully implicit solver implementations.
/// </summary>
namespace FIM.Solver
{
    /// <summary>
    /// The fully implicit solver algorithm.
    /// </summary>
    public class FullyImplicitSolver
    {
        // the entry point

        // this variable is used to store the current time of the simulation run.
        static double currentTime;

        /// <summary>
        /// The entry point of the fully implicit simulation cycle.
        /// </summary>
        /// <param name="data">The <see cref="SimulationData"/>.</param>
        public static void RunSimulation(SimulationData data)
        {
            // size of the model grid.
            int size = data.grid.Length * data.phases.Length;

            // initialize the initial Jacobi matrix.
            double[][] jacobian = new double[size][];
            for (int i = 0; i < jacobian.Length; i++)
            {
                jacobian[i] = new double[size];
            }

            // initalize the empty minus R matrix.
            double[] minusR = new double[size];

            // intialize the empty delta_x matrix. this matrix will store the solution of the Ax=B equation "where A is the Jacobi and B is minus R".
            double[] delta = new double[size];

            // the total simulated time.
            double end_time = data.endingTime;

            // the simulation loop.

            currentTime = 0;
            data.output.Write(currentTime, data, true);

            for (; currentTime < end_time;)
            {
                IterativeSolver(data, jacobian, minusR, delta);
                // this way of updating the loop iterator "currentTime" is used to correctly display current time
                // after the iterations. This is because the time step length may change during the iteration.
                currentTime += data.timeStep;

                //Console.WriteLine(currentTime + ", " + data.grid[0].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.MBE_Gas + ", " + data.wells[0].BHP[1] + ", " + data.wells[0].q_free_gas[0] + ", " + data.wells[0].q_solution_gas[0] + ", " + data.wells[0].q_oil[0]);
                data.output.Write(currentTime, data);
                //Console.ReadKey();
            }

        }

        // the iteration cycle
        private static void IterativeSolver(SimulationData data, double[][] Jacobi, double[] minusR, double[] delta)
        {
            // an array containing previous "index 0" and current "index 1" material balance errors.
            double[] MBE = new double[2];

            // resets time step to its original value at the beginning of each new time step.
            ResetTimeStep(data);

            // reset relaxation factor.
            data.relaxationFactor = data.originalRelaxationFactor;

            // non-linear iteration counter.
            int counter = 0;

            do
            {
                // formulation.
                NumericalPerturbation.CalculateMinusR_Matrix(data, minusR);
                NumericalPerturbation.CalculateJacobi_Matrix(data, minusR, Jacobi);
                WellTerms.Add(data, Jacobi, minusR);

                // solving the equations.
                delta = SolveForDelta(Jacobi, minusR);

                // update properties with better approximations.
                counter = Update(data, MBE, delta, counter);

                // repeat the cycle if the material balance error is larger than the tolerance specified.
                // as long as the repetitions are smaller than the maximum number of non-linear iterations specified.
            } while (MBE[1] >= data.MBE_Tolerance && counter <= data.maximumNonLinearIterations);

            // update properties of the new time step after successful convergence.
            UpdateProperties(data);
        }

        // contains the algorithms of updating properties, convergence checking, relaxation of deltas and time cut off.
        private static int Update(SimulationData data, double[] MBE,double[] deltas, int counter)
        {
            // check convergence errors.
            bool repeat = CheckConvergence(data, MBE);

            // will not repeat the time step. No stagnation or oscillations and the MBE is decreasing.
            if (!repeat && (data.maximumNonLinearIterations - counter) != 1)
            {
                StabilizeNewton(data, deltas, MBE);

                UpdatePropertiesFromDeltas(data, deltas);

                // increment the counter.
                counter += 1;
            }
            else
            {
                // slash the time step.
                data.timeStep *= data.timeStepSlashingFactor;

                // reset relaxation factor.
                data.relaxationFactor = data.originalRelaxationFactor;

                for (int i = 0; i < data.grid.Length; i++)
                {
                    data.grid[i].Reset(data);
                }

                // reset the convergenceError array so that it violates the convergence criteria
                MBE[0] = 0;
                MBE[1] = data.MBE_Tolerance;
                //convergenceError[1] = 0;

                // reset the counter.
                // this way we begin repeating the same time step with all the number of non linear iterations.
                counter = 0;
            }

            return counter;
        }

        // solves the Ax=B equation "where A is the Jacobi and B is minus R" for deltas.
        private static double[] SolveForDelta(double[][] Jacobi, double[] minusR)
        {
            // note that this method while calculating delta_x "from the Jacobi and minus R matrices", it modifies them as a side effect.
            // this should be taken with much care.
            Mathematics.SolveLinearEquation.direct(Jacobi, minusR);
            return minusR;
        }

        // updates the n1 time level properties with better approximation after solving the Ax=B equation.
        private static void UpdatePropertiesFromDeltas(SimulationData data, double[] delta)
        {
            double P, So, Sg, Sw;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                P = data.grid[i].P[1] + delta[counter];
                // assert p is always greater than 0.
                P = P > 0 ? P : 0;

                Sg = data.grid[i].Sg[1] + delta[counter + 1];
                // assert Sg is always greater than 0.
                Sg = Sg > 0 ? Sg : 0;

                Sw = data.grid[i].Sw[1] + delta[counter + 2];
                // assert Sw is always greater than 0.
                Sw = Sw > 0 ? Sw : 0;

                So = 1 - Sw - Sg;

                data.grid[i].UpdateProperties(data, P, Sw, Sg, So, 1);

                counter += data.phases.Length;
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

        // miscellaneous internal helper methods

        // returns the absolute value of the gas material balance error of the whole model.
        private static bool CheckConvergence(SimulationData data, double[] convergenceError)
        {
            bool firstIteration = convergenceError[0] == 0;

            convergenceError[0] = convergenceError[1];

            data.MBE_Oil = Math.Abs(MBE.CheckOil(data));
            data.MBE_Gas = Math.Abs(MBE.CheckGas(data));
            data.MBE_Water = Math.Abs(MBE.CheckWater(data));

            convergenceError[1] = Math.Max(data.MBE_Gas, Math.Max(data.MBE_Oil, data.MBE_Water));

            bool MBE_Increasing = (convergenceError[1] > convergenceError[0]);
            bool slowConvergence = convergenceError[1] / convergenceError[0] > data.maximumMaterialBalanceErrorRatio;

            // if this is not the first iteration and either one of the following conditions is true:
            // 1- material balance error is increasing not decreasing.
            // 2- the rate of convergence is slow.
            if (!firstIteration && (MBE_Increasing || slowConvergence))
            {
                data.relaxationFactor -= data.relaxationFactorDecrement;

                if (data.relaxationFactor < data.minimumRelaxation)
                {
                    return true;
                }

                return false;
            }

            return false;
        }

        /// <summary>
        /// this method detects oscillations or stagnations and applies relaxation if needed.
        /// </summary>
        /// <remarks>
        /// The only time this method triggered was when i ran the model using <see cref="AccumulationTermExpansion"/> near the bubble point pressure.
        /// </remarks>
        /// <param name="data">the <see cref="SimulationData"/>.</param>
        /// <param name="delta_x">The result of solving the Ax=B equation</param>
        /// <param name="convergenceError">
        /// <para>An array containing previous and current iteration material balance error.</para>
        /// <para>This way we can detect oscillations or stagnation based on comparing current and previous iteration values.</para>
        /// </param>
        /// <returns>
        /// <para>True if relaxation fails. This requires repition of the step using a different time step</para>
        /// <para>False if relaxation is applied.</para>
        /// </returns>
        /// <seealso cref="SimulationData.timeStepSlashingFactor"/>.
        private static void StabilizeNewton(SimulationData data, double[] delta_x, double[] convergenceError)
        {

            if (data.relaxationFactor != 1)
            {
                for (int i = 0; i < delta_x.Length; i++)
                {
                    delta_x[i] = data.relaxationFactor * delta_x[i];
                }
            }

        }

        // resets the time step to the original value.
        // this code may be modified to make gradual return from the cut off time step to the original value.
        private static void ResetTimeStep(SimulationData data)
        {
            // this will make sure the time step is reset to the original value
            // it also adjusts the time step sothat the last step will end at the specified ending time.
            data.timeStep = currentTime + data.originalTimeStep > data.endingTime ? data.endingTime - currentTime : data.originalTimeStep;
        }

    }
}
