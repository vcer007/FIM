using System.Collections.Generic;
using System.Linq;
/// <summary>
/// This namespace contains classes used for mathematical calculations.
/// </summary>
namespace FIM.Mathematics
{
    /// <summary>
    /// This class contains method used for solving a set of linear equations.
    /// </summary>
    public class SolveLinearEquation
    {
        /// <summary>
        /// A direct solver that calculates delta x.
        /// </summary>
        /// <remarks>
        /// This method is based on gaussian elimination.
        /// </remarks>
        /// <param name="X">The Jacobi matrix.</param>
        /// <param name="Y">The minus_R array.</param>
        public static void direct(double[][] X, double[] Y)
        {
            int I, J, K, K1, N;
            N = Y.Length;
            for (K = 0; K < N; K++)
            {
                K1 = K + 1;
                for (I = K; I < N; I++)
                {
                    if (X[I][K] != 0)
                    {
                        for (J = K1; J < N; J++)
                        {
                            X[I][J] /= X[I][K];
                        }
                        Y[I] /= X[I][K];
                    }
                }
                for (I = K1; I < N; I++)
                {
                    if (X[I][K] != 0)
                    {
                        for (J = K1; J < N; J++)
                        {
                            X[I][J] -= X[K][J];
                        }
                        Y[I] -= Y[K];
                    }
                }
            }
            for (I = N - 2; I >= 0; I--)
            {
                for (J = N - 1; J >= I + 1; J--)
                {
                    Y[I] -= X[I][J] * Y[J];
                }
            }
        }

        public static void GaussSeidel(FIMSparseMatrix A, double[] X, double[] B)
        {
            double temp = 0;

            int counter;
            int linearIterationsCounter = 0;

            do
            {
                counter = 0;
                linearIterationsCounter++;

                for (int i = 0; i < X.Length; i++)
                {
                    temp = 0;
                    while (counter < A.indicesY.Length && A.indicesY[counter] == i)
                    {
                        temp += A.values[counter] * X[A.indicesX[counter]];
                        counter++;
                    }
                    temp -= A.values[A.indicesDiagonal[i]] * X[i];
                    X[i] = 1 / A.values[A.indicesDiagonal[i]] * (B[i] - temp);
                }
            } while (linearIterationsCounter < 10);

        }

        public static void Orthomin(FIMSparseMatrix A, ref double[] X, double[] B)
        {
            int orthominIterations = 8;
            int size = X.Length;
            double maxError = 0;
            double tolerance = 1E-10;
            int maximumLinearIterations = 1000;
            int counter = 0;
            bool breakLoop = false;

            double a;
            double[] r, temp;

            List<double> b = new List<double>();

            List<double> dotg = new List<double>();
            List<double[]> g = new List<double[]>();
            List<double[]> p = new List<double[]>();

            temp = new double[size];

            r = VectorOperations.Subtract(B, A.Multiply(X));
            p.Add(r.Clone() as double[]);

            do
            {
                for (int i = 0; i < orthominIterations; i++)
                {
                    counter += 1;

                    g.Add(A.Multiply(p[i]));
                    dotg.Add(VectorOperations.Dot(g[i], g[i]));

                    a = VectorOperations.Dot(r, g[i]) / dotg[i];
                    X = VectorOperations.Add(X, VectorOperations.Multiply(a, p[i]));
                    r = VectorOperations.Subtract(r, VectorOperations.Multiply(a, g[i]));
                    maxError = r.Max();
                    breakLoop = maxError <= tolerance;
                    if (breakLoop) break;

                    temp = A.Multiply(r);
                    b.Clear();
                    for (int j = 0; j < g.Count; j++)
                    {
                        b.Add(-VectorOperations.Dot(temp, g[j]) / dotg[j]);
                    }
                    p.Add(VectorOperations.Add(r, VectorOperations.SumMultiply(b, p)));

                }

                if (breakLoop) break;
                //reset
                dotg.Clear();
                g.Clear();
                p.Clear();

                p.Add(r.Clone() as double[]);

            } while (!breakLoop && counter < maximumLinearIterations);

        }
    }
}
