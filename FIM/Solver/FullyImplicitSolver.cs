using FIM.Core;
using FIM.MaterialBalance;
using FIM.Extensions.FullyImplicit;
using FIM.Extensions;

using System;
using FIM.Mathematics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Double;

using System.Linq;
using MathNet.Numerics.Providers.LinearAlgebra.Mkl;

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

        // private members put here globally in the class to reduce the number of GC runs
        static Matrix<double> A;
        static Vector<double> B;
        static Vector<double> X;

        // an array containing previous "index 0" and current "index 1" material balance errors.
        static double[] MBE = new double[2];

        // Create monitor with defined stop criteria
        //static Iterator<double> monitor = new Iterator<double>(iterationCountStopCriterion, residualStopCriterion);

        static IIterativeSolver<double> solver = new GpBiCg();
        static IPreconditioner<double> preconditioner = new MILU0Preconditioner(false);

        /// <summary>
        /// The entry point of the fully implicit simulation cycle.
        /// </summary>
        /// <param name="data">The <see cref="SimulationData"/>.</param>
        public static void RunSimulation(SimulationData data)
        {
            // size of the model grid.
            int size = data.grid.Length * data.phases.Length;

            var jacobi = new MathNet.Numerics.LinearAlgebra.Double.SparseMatrix(size, size);

            //MathNetSparseMatrixInitializer.Initialize(data);

            //jacobi.InitializeStructure(data, 1);

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
                var time_stamp = DateTime.Now;

                int counter = IterativeSolver(data, jacobi, minusR, delta);

                var duration = DateTime.Now - time_stamp;
                // this way of updating the loop iterator "currentTime" is used to correctly display current time
                // after the iterations. This is because the time step length may change during the iteration.
                currentTime += data.timeStep;

                //Console.WriteLine(currentTime + ", " + data.grid[0].P[0] + ", " + data.grid[0].Sg[0] + ", " + data.MBE_Gas + ", " + data.wells[0].BHP[1] + ", " + data.wells[0].q_free_gas[0] + ", " + data.wells[0].q_solution_gas[0] + ", " + data.wells[0].q_oil[0]);
                data.output.Write(currentTime, data, false, counter, (int)duration.TotalMilliseconds);
                //Console.ReadKey();
            }

        }

        // the iteration cycle
        private static int IterativeSolver(SimulationData data, MathNet.Numerics.LinearAlgebra.Double.SparseMatrix Jacobi, double[] minusR, double[] delta)
        {
            //use the native Intel MKL solver "which is several times faster"

            Control.UseManaged();
            //Control.UseNativeMKL();
            //Control.NativeProviderPath = @"C:\MKL";

            MBE[0] = 0; MBE[1] = 0;

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
                //var test = Jacobi.Storage.Enumerate().ToList();
                //var test2 = test.Where(x => x != 0);
                //var test3 = test2.Where(x => x > -1E-7 && x < 1E-7);
                //var test4 = test3.Select(x => Math.Abs(x)).Min();
                WellTerms.Add(data, Jacobi, minusR);

                // solving the equations.
                //SparseMatrix matrix = SparseMatrix.FromJagged(Jacobi);
                //for (int i = 0; i < delta.Length; i++)
                //{
                //    delta[i] = Global.EPSILON;
                //}
                //Mathematics.SolveLinearEquation.Orthomin(matrix, ref delta, minusR);

                //A = Matrix<double>.Build.SparseOfRowArrays(Jacobi);
                B = Vector<double>.Build.Dense(minusR);
                X = Vector<double>.Build.Dense(delta);
                //var X = Vector<double>.Build.DenseOfArray(delta);

                // Stop calculation if 1000 iterations reached during calculation
                IIterationStopCriterion<double> iterationCountStopCriterion = new IterationCountStopCriterion<double>(1000);
                // Stop calculation if residuals are below 1E-10 --> the calculation is considered converged
                IIterationStopCriterion<double> residualStopCriterion = new ResidualStopCriterion<double>(1e-5);

                //delta = A.SolveIterative(B, solver, preconditioner, iterationCountStopCriterion, residualStopCriterion).AsArray();
                var iterativeSolverStatus = Jacobi.TrySolveIterative(B, X, solver, preconditioner, iterationCountStopCriterion, residualStopCriterion);
                delta = X.AsArray();
                //delta = Jacobi.Solve(B).ToArray();

                ////delta = (A.Inverse() * B).ToArray();
                //delta = A.Solve(B).ToArray();
                //delta = A.SolveIterative(B, solver, monitor).ToArray();
                //delta = A.LU().Solve(B).ToArray();

                //delta = SolveForDelta(Jacobi, minusR);

                // update properties with better approximations.
                counter = Update(data, MBE, delta, counter);
                //if (iterativeSolverStatus == IterationStatus.Converged)
                //{
                //}
                //else
                //{
                //    for (int i = 0; i < data.grid.Length; i++)
                //    {
                //        data.grid[i].Reset(data);
                //    }
                //    data.timeStep *= data.timeStepSlashingFactor;
                //    counter = 0;
                //    MBE[0] = 0;
                //    MBE[1] = data.MBE_Tolerance;
                //}

                //Console.WriteLine(MBE[1]);

                // repeat the cycle if the material balance error is larger than the tolerance specified.
                // as long as the repetitions are smaller than the maximum number of non-linear iterations specified.
            } while (MBE[1] >= data.MBE_Tolerance && counter <= data.maximumNonLinearIterations);

            // update properties of the new time step after successful convergence.
            UpdateProperties(data);

            return counter;

        }

        // contains the algorithms of updating properties, convergence checking, relaxation of deltas and time cut off.
        private static int Update(SimulationData data, double[] MBE,double[] deltas, int counter)
        {
            //StabilizeNewton(data, deltas, MBE);
            UpdatePropertiesFromDeltas(data, deltas);

            // check convergence errors.
            bool repeat = CheckConvergence(data, MBE);

            // will not repeat the time step. No stagnation or oscillations and the MBE is decreasing.
            if (!repeat && (data.maximumNonLinearIterations - counter) != 1)
            {

                //UpdatePropertiesFromDeltas(data, deltas);

                // increment the counter.
                counter += 1;
            }
            else
            {
                // slash the time step.
                data.timeStep *= data.timeStepSlashingFactor;

                //Console.WriteLine($"The new time-step is {data.timeStep}");

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
            double temp;
            for (int i = 0; i < data.grid.Length; i++)
            {
                temp = data.grid[i].P[1] + delta[counter];
                P = temp > 0 ? temp : data.grid[i].P[1];
                // assert p is always greater than 0.
                //P = P > 0 ? P : 0;

                temp = data.grid[i].Sg[1] + delta[counter + 1];
                Sg = temp > 0 && temp < 1 ? temp : data.grid[i].Sg[1];
                Sg = Sg > 1 ? 1 : Sg;
                //Sg = temp > 1 ? 1 : (temp < 0 ? 0 : temp);

                temp = data.grid[i].Sw[1] + delta[counter + 2];
                Sw = temp > 0 && temp < 1 ? temp : data.grid[i].Sw[1];
                Sw = Sw > 1 ? 1 : Sw;
                //Sw = Sw < data.pvt.connateWaterSaturation ? data.pvt.connateWaterSaturation : Sw;
                //Sw = temp > 1 ? 1 : (temp < 0 ? 0 : temp);

                temp = 1 - Sw - Sg;
                So = temp > 0 && temp < 1 ? temp : data.grid[i].So[1];

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

            data.MBE_Oil = Math.Abs(MaterialBalance.MBE.CheckOil(data));
            data.MBE_Gas = Math.Abs(MaterialBalance.MBE.CheckGas(data));
            data.MBE_Water = Math.Abs(MaterialBalance.MBE.CheckWater(data));

            convergenceError[1] = Math.Max(data.MBE_Gas, Math.Max(data.MBE_Oil, data.MBE_Water));

            bool MBE_Increasing = (convergenceError[1] > convergenceError[0]);
            bool slowConvergence = convergenceError[1] / convergenceError[0] > data.maximumMaterialBalanceErrorRatio;


            if (!firstIteration && MBE_Increasing)
            {
                return true;
            }

            if (!firstIteration && slowConvergence)
            {
                data.relaxationFactor -= data.relaxationFactorDecrement;

                if (data.relaxationFactor < data.minimumRelaxation)
                {
                    return true;
                }
            }

            return false;

            //// if this is not the first iteration and either one of the following conditions is true:
            //// 1- material balance error is increasing not decreasing.
            //// 2- the rate of convergence is slow.
            //if (!firstIteration && (MBE_Increasing || slowConvergence))
            //{
            //    data.relaxationFactor -= data.relaxationFactorDecrement;

            //    //Console.WriteLine($"The new relaxation factor is {data.relaxationFactor}");

            //    if (data.relaxationFactor < data.minimumRelaxation)
            //    {
            //        return true;
            //    }

            //    return false;
            //}

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
