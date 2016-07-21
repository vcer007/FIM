using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Extensions
{
    class WellTerms
    {
        public static void add(SimulationData data, double[][] Jacobi_Matrix, double[] minusR_Matrix)
        {
            int index = 0;

            for (int i = 0; i < data.wells.Length; i++)
            {
                index = 3 * data.wells[i].index;

                minusR_Matrix[index] += data.wells[i].q_oil[1];
                minusR_Matrix[index + 1] += data.wells[i].q_free_gas[1] + data.wells[i].q_solution_gas[1];
                minusR_Matrix[index + 2] += data.wells[i].q_water[1];
            }

            for (int i = 0; i < data.wells.Length; i++)
            {
                index = 3 * data.wells[i].index;

                Jacobi_Matrix[index][index] -= data.wells[i].dq_oil_dP;
                Jacobi_Matrix[index][index + 1] -= data.wells[i].dq_oil_dSg;
                Jacobi_Matrix[index][index + 2] -= data.wells[i].dq_oil_dSw;

                Jacobi_Matrix[index + 1][index] -= (data.wells[i].dq_free_gas_dP + data.wells[i].dq_solution_gas_dP);
                Jacobi_Matrix[index + 1][index + 1] -= (data.wells[i].dq_free_gas_dSg + data.wells[i].dq_solution_gas_dSg);
                Jacobi_Matrix[index + 1][index + 2] -= (data.wells[i].dq_free_gas_dSw + data.wells[i].dq_solution_gas_dSw);

                Jacobi_Matrix[index + 2][index] -= data.wells[i].dq_water_dP;
                Jacobi_Matrix[index + 2][index + 1] -= data.wells[i].dq_water_dSg;
                Jacobi_Matrix[index + 2][index + 2] -= data.wells[i].dq_water_dSw;
            }
        }
    }
}
