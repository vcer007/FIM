using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FIM.Core;

namespace FIM.Mathematics
{
    class SolveLinearEquation
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
    }
}
