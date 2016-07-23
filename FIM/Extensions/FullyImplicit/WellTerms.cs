using FIM.Core;

namespace FIM.Extensions.FullyImplicit
{
    /// <summary>
    /// This class handles the addition of the well terms to the Jacobi and the minus R matrices.
    /// </summary>
    public static class WellTerms
    {
        /// <summary>
        /// Add the well terms "flow rates and derivatives" of each phase to their respective position in the Jacobi and the minus R matrices.
        /// </summary>
        /// <param name="data"></param>
        /// <param name="Jacobi"></param>
        /// <param name="minusR"></param>
        public static void add(SimulationData data, double[][] Jacobi, double[] minusR)
        {
            int index = 0;

            for (int i = 0; i < data.wells.Length; i++)
            {
                index = 3 * data.wells[i].index;

                minusR[index] += data.wells[i].q_oil[1];
                minusR[index + 1] += data.wells[i].q_free_gas[1] + data.wells[i].q_solution_gas[1];
                minusR[index + 2] += data.wells[i].q_water[1];
            }

            for (int i = 0; i < data.wells.Length; i++)
            {
                index = 3 * data.wells[i].index;

                Jacobi[index][index] -= data.wells[i].dq_oil_dP;
                Jacobi[index][index + 1] -= data.wells[i].dq_oil_dSg;
                Jacobi[index][index + 2] -= data.wells[i].dq_oil_dSw;

                Jacobi[index + 1][index] -= (data.wells[i].dq_free_gas_dP + data.wells[i].dq_solution_gas_dP);
                Jacobi[index + 1][index + 1] -= (data.wells[i].dq_free_gas_dSg + data.wells[i].dq_solution_gas_dSg);
                Jacobi[index + 1][index + 2] -= (data.wells[i].dq_free_gas_dSw + data.wells[i].dq_solution_gas_dSw);

                Jacobi[index + 2][index] -= data.wells[i].dq_water_dP;
                Jacobi[index + 2][index + 1] -= data.wells[i].dq_water_dSg;
                Jacobi[index + 2][index + 2] -= data.wells[i].dq_water_dSw;
            }
        }
    }
}
