using FIM.Core;
using FIM.Well;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Extensions.IMPES
{
    /// <summary>
    /// This class contains helper method for formulating the equations used to solve the IMPES problem.
    /// </summary>
    /// <remarks>
    /// Based on equation 35, page 11, chapter 4, Wattenberger PETE 603 notes.
    /// </remarks>
    class IMPES_Preparation
    {

        public static void CalculateConstantsMatrix(SimulationData data, double[] constants)
        {
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                constants[i] = -block.Vp[0] / (Global.a * data.timeStep) * block.GetCompressibilityTotal(data) * block.P[0];
            }
        }

        public static void CalculatePressureCoefficientsMatrix(SimulationData data, double[][] pressureCoefficients)
        {
            // set the matrix elements to zero.
            // this is necessary as the direct solver used alters the values of the matrix as a side-effect.
            for (int i = 0; i < pressureCoefficients.Length; i++)
            {
                for (int j = 0; j < pressureCoefficients[i].Length; j++)
                {
                    pressureCoefficients[i][j] = 0;
                }
            }

            BaseBlock block, neighbor_block, upstream_block, downstream_block;
            double transmissibility, kr, viscosity, B, RsO;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];


                // this block
                pressureCoefficients[i][block.index] += -block.Vp[0] / (Global.a * data.timeStep) *  block.GetCompressibilityTotal(data);

                // neighboring blocks

                // oil phase
                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[j]];
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

                    pressureCoefficients[i][block.index] += -1 * block.Bo[1] * transmissibility * kr / (viscosity * B) /** block.P[1]*/;

                    pressureCoefficients[i][neighbor_block.index] += block.Bo[1] * transmissibility * kr / (viscosity * B) /** neighbor_block.P[1]*/;
                }

                // free gas
                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[j]];
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

                    pressureCoefficients[i][block.index] += -1 * block.Bg[1] * transmissibility * kr / (viscosity * B) /** block.P[1]*/;

                    pressureCoefficients[i][neighbor_block.index] += block.Bg[1] * transmissibility * kr / (viscosity * B) /** neighbor_block.P[1]*/;
                }

                // solution gas
                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[j]];
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

                    pressureCoefficients[i][block.index] += -1 * block.Bg[1] * RsO * transmissibility * kr / (viscosity * B) /** block.P[1]*/;

                    pressureCoefficients[i][neighbor_block.index] += block.Bg[1] * RsO * transmissibility * kr / (viscosity * B) /** neighbor_block.P[1]*/;

                    pressureCoefficients[i][block.index] += block.Rso[1] * block.Bg[1] * transmissibility * kr / (viscosity * B) /** block.P[1]*/;

                    pressureCoefficients[i][neighbor_block.index] += -1 * block.Rso[1] * block.Bg[1] * transmissibility * kr / (viscosity * B) /** neighbor_block.P[1]*/;
                }

                // water phase
                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[j]];
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

                    pressureCoefficients[i][block.index] += -1 * block.Bw[1] * transmissibility * kr / (viscosity * B) /** block.P[1]*/;

                    pressureCoefficients[i][neighbor_block.index] += block.Bw[1] * transmissibility * kr / (viscosity * B) /** neighbor_block.P[1]*/;
                }
            }
        }

        public static void AddWellTerms(SimulationData data, double[][] pressureCoefficients, double[] constants)
        {
            int index;
            double temp;
            BaseBlock block;
            Well.BaseWell well;

            //double q_oil, q_water, q_gas, drawDown;

            for (int i = 0; i < data.wells.Length; i++)
            {
                well = data.wells[i];
                index = well.index;
                block = data.grid[index];

                
                if (well.control == Global.WellControl.OilRate)
                {
                    constants[index] += well.q_oil[0] * block.Bo[1] + well.q_water[0] * block.Bw[1] + well.q_free_gas[0] * block.Bg[1];
                }
                else if (well.control == Global.WellControl.BHP)
                {
                    if (well.method == Global.WellRateCalculation.Explicit)
                    {
                        constants[index] += well.q_oil[0] * block.Bo[1] + well.q_water[0] * block.Bw[1] + well.q_free_gas[0] * block.Bg[1];
                    }
                    else
                    {
                        // oil
                        temp = BaseWell.GetWellMobility(block.Kro[0], block.viscosityOil[0], block.Bo[0]);
                        pressureCoefficients[index][index] -= block.Bo[1] * temp;
                        constants[index] -= block.Bo[1] * temp * well.BHP[1];

                        // solution gas
                        //LHS_Matrix[index][index] -= block.Bg[1] * block.Rso[1] * temp * block.P[1];
                        //RHS_Matrix[index] -= block.Bg[1] * block.Rso[1] * temp * well.specifiedMinimumBHP;

                        // water
                        temp = BaseWell.GetWellMobility(block.Krw[0], block.viscosityWater[0], block.Bw[0]);
                        pressureCoefficients[index][index] -= block.Bw[1] * temp;
                        constants[index] -= block.Bw[1] * temp * well.BHP[1];

                        // gas
                        temp = BaseWell.GetWellMobility(block.Krg[0], block.viscosityGas[0], block.Bg[0]);
                        pressureCoefficients[index][index] -= block.Bg[1] * temp;
                        constants[index] -= block.Bg[1] * temp * well.BHP[1];
                    }
                    
                }

            }
        }

    }
}
