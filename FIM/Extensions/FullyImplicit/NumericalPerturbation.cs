using FIM.Core;
using MathNet.Numerics.LinearAlgebra.Double;

using FIM.Mathematics;
using MathNet.Numerics.LinearAlgebra.Storage;
using System.Linq;
using System.Threading.Tasks;

/// <summary>
/// The name space that organizes extension methods to <see cref="BaseBlock"/> specific to fully implicit simulation.
/// </summary>
namespace FIM.Extensions.FullyImplicit
{
    /// <summary>
    /// This class contains extension methods to the <see cref="BaseBlock"/> to add support for perturbation.
    /// </summary>
    public static class NumericalPerturbation
    {

        // calculate the value of the residual equation for any phase.
        private static double CalculateR(this BaseBlock block, SimulationData data, Global.Phase phase)
        {
            bool GravityCalculation = data.Gravity;

            // Note : Kro, Krw, Krg are all calculated based on Sg.

            double gravity = 0;

            BaseBlock upstream_block, downstream_block, neighbor_block;
            double transmissibility;

            double Kr, B, viscosity, Rso, Rvo;
            double R = 0;

            double temp, accumulation_term;


            if (phase == Global.Phase.Oil)
            {
                for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                {
                    if (block.neighborBlocksIndices[i] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[i]];
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (GravityCalculation)
                    {
                        gravity = data.pvt.GetAverageOilGravity(block, neighbor_block);
                    }

                    transmissibility = block.transmissibility_list[i];

                    if (GetOilPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Kro[1];
                    B = 0.5 * (block.Bo[1] + neighbor_block.Bo[1]);
                    viscosity = 0.5 * (block.viscosityOil[1] + neighbor_block.viscosityOil[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1] - delta_z * gravity);
                    block.transmissibility_terms_oil[i] = temp;

                    if (data.vaporizedOilPresent)
                    {
                        Kr = upstream_block.Krg[1];
                        B = 0.5 * (block.Bg[1] + neighbor_block.Bg[1]);
                        viscosity = 0.5 * (block.viscosityGas[1] + neighbor_block.viscosityGas[1]);
                        Rvo = 0.5 * (block.Rvo[1] + neighbor_block.Rvo[1]);

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageGasGravity(block, neighbor_block);
                        }

                        temp += Rvo * transmissibility * Kr / (viscosity * B) * (neighbor_block.GetPg(1, 1) - block.GetPg(1, 1) - delta_z * gravity);
                    }

                    R += temp;

                }

                accumulation_term = 1 / (Global.a * data.timeStep) * ((block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));

                if (data.vaporizedOilPresent)
                {
                    accumulation_term += 1 / (Global.a * data.timeStep) * ((block.Rvo[1] * block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Rvo[0] * block.Vp[0] * block.Sg[0] / block.Bg[0]));
                }

                block.accumulation_term_oil = accumulation_term;

                R -= accumulation_term;
            }
            else if (phase == Global.Phase.Gas)
            {
                for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                {
                    if (block.neighborBlocksIndices[i] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[i]];
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (GravityCalculation)
                    {
                        gravity = data.pvt.GetAverageGasGravity(block, neighbor_block);
                    }

                    transmissibility = block.transmissibility_list[i];

                    if (GetGasPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Krg[1];
                    B = 0.5 * (block.Bg[1] + neighbor_block.Bg[1]);
                    viscosity = 0.5 * (block.viscosityGas[1] + neighbor_block.viscosityGas[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbor_block.GetPg(1, 1) - block.GetPg(1, 1) - delta_z * gravity);

                    if (data.solubleGasPresent)
                    {
                        B = 0.5 * (block.Bo[1] + neighbor_block.Bo[1]);
                        viscosity = 0.5 * (block.viscosityOil[1] + neighbor_block.viscosityOil[1]);
                        Rso = 0.5 * (block.Rso[1] + neighbor_block.Rso[1]);

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageOilGravity(block, neighbor_block);
                        }

                        if (GetOilPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                        {
                            upstream_block = block;
                            downstream_block = neighbor_block;
                        }
                        else
                        {
                            upstream_block = neighbor_block;
                            downstream_block = block;
                        }

                        Kr = upstream_block.Kro[1];

                        temp += Rso * transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1] - delta_z * gravity);
                    }

                    block.transmissibility_terms_gas[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / (Global.a * data.timeStep) * ((block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                // check for presence of soluble_gas in simulation_data
                if (data.solubleGasPresent)
                {
                    accumulation_term += 1 / (Global.a * data.timeStep) * ((block.Rso[1] * block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                block.accumulation_term_gas = accumulation_term;

                R -= accumulation_term;
            }
            else// if (phase == Global.Phase.Water)
            {
                for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                {
                    if (block.neighborBlocksIndices[i] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighborBlocksIndices[i]];

                    if (GravityCalculation)
                    {
                        gravity = data.pvt.GetAverageWaterGravity(block, neighbor_block);
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

                    transmissibility = block.transmissibility_list[i];

                    if (GetWaterPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Krw[1];
                    B = 0.5 * (block.Bw[1] + neighbor_block.Bw[1]);
                    viscosity = 0.5 * (block.viscosityWater[1] + neighbor_block.viscosityWater[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1] - delta_z * gravity);
                    block.transmissibility_terms_water[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / (Global.a * data.timeStep) * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                R -= accumulation_term;
            }
            
            return R;
        }

        // perturb the residual equation "calculate at independent variable + epsilon".
        // the perturbed value will then be subtracted from the original R value and divided by epsilon to get the derivative numerically.
        private static double Perturb(this BaseBlock block, SimulationData data, Global.Phase equation_phase, int neighbour_block_index, Global.Variable variable)
        {
            double gravity = 0;

            double R_plus = 0;

            double transmissibility_temp = 0;

            double transmissibility_term = 0;
            double accumulation_term = 0;

            double kr = 1, viscosity, B, Rso, Rvo;
            BaseBlock upstream_block, downstream_block, neighbor_block;

            double transmissiblity;

            #region Variable Dependencies
            int pd = 1, npd = 1, swd = 1, nswd = 1, sgd = 1, nsgd = 1;
            bool this_block = true;

            if (neighbour_block_index == -1)
            {
                this_block = true;
            }
            else
            {
                this_block = false;
            }

            if (variable == Global.Variable.Pressure)
            {
                if (this_block)
                {
                    pd = 2;
                }
                else
                {
                    npd = 2;
                }
            }
            else if (variable == Global.Variable.SaturationGas)
            {
                if (this_block)
                {
                    sgd = 2;
                }
                else
                {
                    nsgd = 2;
                }
            }
            else if (variable == Global.Variable.SaturationWater)
            {
                if (this_block)
                {
                    swd = 2;
                }
                else
                {
                    nswd = 2;
                }
            }
            #endregion

            bool GravityCalculation = data.Gravity;

            if (equation_phase == Global.Phase.Oil)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                    {
                        if (block.neighborBlocksIndices[i] == -1)
                        {
                            continue;
                        }

                        neighbor_block = data.grid[block.neighborBlocksIndices[i]];
                        double delta_z = neighbor_block.Depth - block.Depth;

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageOilGravity(block, neighbor_block);
                        }

                        kr = GetKr(Global.Phase.Oil, block, neighbor_block, variable, this_block, delta_z, gravity, data);

                        B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd] - delta_z * gravity);

                        if (data.vaporizedOilPresent)
                        {

                            Rvo = 0.5 * (block.Rvo[pd] + neighbor_block.Rvo[npd]);
                            B = 0.5 * (block.Bg[pd] + neighbor_block.Bg[npd]);
                            viscosity = 0.5 * (block.viscosityGas[pd] + neighbor_block.viscosityGas[npd]);

                            if (GravityCalculation)
                            {
                                gravity = data.pvt.GetAverageGasGravity(block, neighbor_block);
                            }

                            if (neighbor_block.GetPg(1, 1) - block.GetPg(1, 1) - delta_z * gravity <= 0)
                            {
                                upstream_block = block;
                                downstream_block = neighbor_block;

                                kr = upstream_block.Krg[sgd];
                            }
                            else
                            {
                                upstream_block = neighbor_block;
                                downstream_block = block;

                                kr = upstream_block.Krg[nsgd];
                            }

                            transmissibility_temp += Rvo * transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPg(npd, nsgd) - block.GetPg(pd, sgd) - delta_z * gravity);
                        }

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.timeStep * Global.a) * ((block.Vp[pd] * (1 - block.Sg[sgd] - block.Sw[swd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0]));

                    if (data.vaporizedOilPresent)
                    {
                        accumulation_term += 1 / (data.timeStep * Global.a) * (block.Rvo[pd] * (block.Vp[pd] * block.Sg[sgd] / block.Bg[pd]) - block.Rvo[0] * (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                    }
                }
                else
                {
                    // transmissibility term
                    neighbor_block = data.grid[block.neighborBlocksIndices[neighbour_block_index]];
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (GravityCalculation)
                    {
                        gravity = data.pvt.GetAverageOilGravity(block, neighbor_block);
                    }

                    kr = GetKr(Global.Phase.Oil, block, neighbor_block, variable, this_block, delta_z, gravity, data);

                    B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                    viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd] - delta_z * gravity);

                    if (data.vaporizedOilPresent)
                    {

                        Rvo = 0.5 * (block.Rvo[pd] + neighbor_block.Rvo[npd]);
                        B = 0.5 * (block.Bg[pd] + neighbor_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosityGas[pd] + neighbor_block.viscosityGas[npd]);

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageGasGravity(block, neighbor_block);
                        }

                        if (neighbor_block.GetPg(1, 1) - block.GetPg(1, 1) - delta_z * gravity <= 0)
                        {
                            upstream_block = block;
                            downstream_block = neighbor_block;

                            kr = upstream_block.Krg[sgd];
                        }
                        else
                        {
                            upstream_block = neighbor_block;
                            downstream_block = block;

                            kr = upstream_block.Krg[nsgd];
                        }

                        transmissibility_temp += Rvo * transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPg(npd, nsgd) - block.GetPg(pd, sgd) - delta_z * gravity);
                    }

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += block.transmissibility_terms_oil[i];
                    }
                    // accumulation term
                    accumulation_term = block.accumulation_term_oil;

                }
            }
            else if (equation_phase == Global.Phase.Gas)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                    {
                        if (block.neighborBlocksIndices[i] == -1)
                        {
                            continue;
                        }

                        neighbor_block = data.grid[block.neighborBlocksIndices[i]];
                        double delta_z = neighbor_block.Depth - block.Depth;

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageGasGravity(block, neighbor_block);
                        }

                        if (GetGasPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                        {
                            upstream_block = block;
                            downstream_block = neighbor_block;

                            kr = upstream_block.Krg[sgd];
                        }
                        else
                        {
                            upstream_block = neighbor_block;
                            downstream_block = block;

                            kr = upstream_block.Krg[nsgd];
                        }

                        B = 0.5 * (block.Bg[pd] + neighbor_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosityGas[pd] + neighbor_block.viscosityGas[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPg(npd, nsgd) - block.GetPg(pd, sgd) - delta_z * gravity);

                        if (data.solubleGasPresent)
                        {

                            Rso = 0.5 * (block.Rso[pd] + neighbor_block.Rso[npd]);
                            B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                            viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);

                            if (GravityCalculation)
                            {
                                gravity = data.pvt.GetAverageOilGravity(block, neighbor_block);
                            }

                            kr = GetKr(Global.Phase.Oil, block, neighbor_block, variable, this_block, delta_z, gravity, data);

                            transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd] - delta_z * gravity);
                        }


                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.timeStep * Global.a) * ((block.Vp[pd] * block.Sg[sgd] / block.Bg[pd]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));

                    if (data.solubleGasPresent)
                    {
                        accumulation_term += 1 / (data.timeStep * Global.a) * (block.Rso[pd] * (block.Vp[pd] * (1 - block.Sg[sgd] - block.Sw[swd]) / block.Bo[pd]) - block.Rso[0] * (block.Vp[0] * block.So[0] / block.Bo[0]));
                    }

                }
                else
                {
                    // transmissibility term
                    neighbor_block = data.grid[block.neighborBlocksIndices[neighbour_block_index]];
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (GravityCalculation)
                    {
                        gravity = data.pvt.GetAverageGasGravity(block, neighbor_block);
                    }

                    if (GetGasPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;

                        kr = upstream_block.Krg[sgd];
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;

                        kr = upstream_block.Krg[nsgd];
                    }

                    B = 0.5 * (block.Bg[pd] + neighbor_block.Bg[npd]);
                    viscosity = 0.5 * (block.viscosityGas[pd] + neighbor_block.viscosityGas[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPg(npd, nsgd) - block.GetPg(pd, sgd) - delta_z * gravity);

                    if (data.solubleGasPresent)
                    {

                        Rso = 0.5 * (block.Rso[pd] + neighbor_block.Rso[npd]);
                        B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageOilGravity(block, neighbor_block);
                        }

                        kr = GetKr(Global.Phase.Oil, block, neighbor_block, variable, this_block, delta_z, gravity, data);

                        transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd] - delta_z * gravity);
                    }

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += block.transmissibility_terms_gas[i];
                    }
                    // accumulation term
                    accumulation_term = block.accumulation_term_gas;
                }
            }
            else if (equation_phase == Global.Phase.Water)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                    {
                        if (block.neighborBlocksIndices[i] == -1)
                        {
                            continue;
                        }

                        neighbor_block = data.grid[block.neighborBlocksIndices[i]];
                        double delta_z = neighbor_block.Depth - block.Depth;

                        if (GravityCalculation)
                        {
                            gravity = data.pvt.GetAverageWaterGravity(block, neighbor_block);
                        }

                        if (GetWaterPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                        {
                            upstream_block = block;
                            downstream_block = neighbor_block;

                            kr = upstream_block.Krw[swd];
                        }
                        else
                        {
                            upstream_block = neighbor_block;
                            downstream_block = block;

                            kr = upstream_block.Krw[nswd];
                        }

                        B = 0.5 * (block.Bw[pd] + neighbor_block.Bw[npd]);
                        viscosity = 0.5 * (block.viscosityWater[pd] + neighbor_block.viscosityWater[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPw(npd, nswd) - block.GetPw(pd, swd) - delta_z * gravity);

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.timeStep * Global.a) * ((block.Vp[pd] * block.Sw[swd] / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));

                }
                else
                {
                    // transmissibility term
                    neighbor_block = data.grid[block.neighborBlocksIndices[neighbour_block_index]];
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (GravityCalculation)
                    {
                        gravity = data.pvt.GetAverageWaterGravity(block, neighbor_block);
                    }

                    if (GetWaterPotentialDifference(neighbor_block, block, delta_z, gravity) <= 0)
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;

                        kr = upstream_block.Krw[swd];
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;

                        kr = upstream_block.Krw[nswd];
                    }

                    B = 0.5 * (block.Bw[pd] + neighbor_block.Bw[npd]);
                    viscosity = 0.5 * (block.viscosityWater[pd] + neighbor_block.viscosityWater[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPw(npd, nswd) - block.GetPw(pd, swd) - delta_z * gravity);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighborBlocksIndices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += block.transmissibility_terms_water[i];
                    }
                    // accumulation term
                    accumulation_term = block.accumulation_term_water;
                }
            }


            R_plus = transmissibility_term - accumulation_term;
            return R_plus;
        }

        /// <summary>
        /// Generates the minus R column array.
        /// </summary>
        /// <param name="data"></param>
        /// <param name="minusR"></param>
        /// <seealso cref="Solver.FullyImplicitSolver"/>
        /// <seealso cref="WellTerms"/>
        public static void CalculateMinusR_Matrix(SimulationData data, double[] minusR)
        {

            BaseBlock block;
            int counter = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                counter = data.phases.Length * i;

                minusR[counter] = -1 * block.CalculateR(data, Global.Phase.Oil);
                minusR[counter + 1] = -1 * block.CalculateR(data, Global.Phase.Gas);
                minusR[counter + 2] = -1 * block.CalculateR(data, Global.Phase.Water);
            }
        }

        /// <summary>
        /// Generates the Jacobi matrix.
        /// </summary>
        /// <param name="data"></param>
        /// <param name="minusR"></param>
        /// <param name="Jacobi"></param>
        /// <seealso cref="Solver.FullyImplicitSolver"/>
        /// <seealso cref="WellTerms"/>
        public static void CalculateJacobi_Matrix(SimulationData data, double[] minusR, SparseMatrix Jacobi)
        {
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                int counter = data.phases.Length * i;

                #region Oil

                for (int j = -1; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (j == -1 || block.neighborBlocksIndices[j] >= 0)
                    {
                        int place = j >= 0 ? data.phases.Length * block.neighborBlocksIndices[j] : data.phases.Length * block.index;
                        // with respect to P
                        var temp = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.Pressure) + minusR[counter]) / Global.EPSILON_P;
                        Jacobi[counter, place] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sg
                        double epsilon_sg = Global.EPSILON_S;
                        if (j == -1)
                        {
                            epsilon_sg = block.epsilon_sg;
                        }
                        else
                        {
                            epsilon_sg = data.grid[block.neighborBlocksIndices[j]].epsilon_sg;
                        }
                        temp = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.SaturationGas) + minusR[counter]) / epsilon_sg;
                        Jacobi[counter, place + 1] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sw
                        temp = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.SaturationWater) + minusR[counter]) / Global.EPSILON_S;
                        Jacobi[counter, place + 2] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;
                    }
                }
                #endregion
                #region Gas

                for (int j = -1; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (j == -1 || block.neighborBlocksIndices[j] >= 0)
                    {
                        int place = j >= 0 ? data.phases.Length * block.neighborBlocksIndices[j] : data.phases.Length * block.index;
                        // with respect to P
                        var temp = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.Pressure) + minusR[counter + 1]) / Global.EPSILON_P;
                        Jacobi[counter + 1, place] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sg
                        double epsilon_sg = Global.EPSILON_S;
                        if (j == -1)
                        {
                            epsilon_sg = block.epsilon_sg;
                        }
                        else
                        {
                            epsilon_sg = data.grid[block.neighborBlocksIndices[j]].epsilon_sg;
                        }
                        temp = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.SaturationGas) + minusR[counter + 1]) / epsilon_sg;
                        Jacobi[counter + 1, place + 1] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sw
                        temp = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.SaturationWater) + minusR[counter + 1]) / Global.EPSILON_S;
                        Jacobi[counter + 1, place + 2] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;
                    }
                }
                #endregion
                #region Water

                for (int j = -1; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (j == -1 || block.neighborBlocksIndices[j] >= 0)
                    {
                        int place = j >= 0 ? data.phases.Length * block.neighborBlocksIndices[j] : data.phases.Length * block.index;
                        // with respect to P
                        var temp = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.Pressure) + minusR[counter + 2]) / Global.EPSILON_P;
                        Jacobi[counter + 2, place] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sg
                        double epsilon_sg = Global.EPSILON_S;
                        if (j == -1)
                        {
                            epsilon_sg = block.epsilon_sg;
                        }
                        else
                        {
                            epsilon_sg = data.grid[block.neighborBlocksIndices[j]].epsilon_sg;
                        }
                        temp = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.SaturationGas) + minusR[counter + 2]) / epsilon_sg;
                        Jacobi[counter + 2, place + 1] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sw
                        temp = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.SaturationWater) + minusR[counter + 2]) / Global.EPSILON_S;
                        Jacobi[counter + 2, place + 2] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;
                    }
                }
                #endregion
            }
        }

        public static double GetKr(Global.Phase phase, BaseBlock block, BaseBlock neighbor, Global.Variable diff_variable, bool thisBlock, double delta_z, double gravity, SimulationData data)
        {
            double potential_difference = 0;

            switch (phase)
            {
                case Global.Phase.Water:
                    potential_difference = GetWaterPotentialDifference(neighbor, block, delta_z, gravity);
                    break;
                case Global.Phase.Oil:
                    potential_difference = GetOilPotentialDifference(neighbor, block, delta_z, gravity);
                    break;
                case Global.Phase.Gas:
                    potential_difference = GetGasPotentialDifference(neighbor, block, delta_z, gravity);
                    break;
                default:
                    break;
            }


            BaseBlock upstream, downstream;
            if (potential_difference <= 0)
            {
                upstream = block;
                downstream = neighbor;
            }
            else
            {
                upstream = neighbor;
                downstream = block;
            }

            double Kr;
            switch (phase)
            {
                case Global.Phase.Water:
                    if (thisBlock && upstream == block && diff_variable == Global.Variable.SaturationWater
                        ||
                        !thisBlock && upstream == neighbor && diff_variable == Global.Variable.SaturationWater)
                    {
                        Kr = upstream.Krw[2];
                    }
                    else
                    {
                        Kr = upstream.Krw[1];
                    }
                    break;
                case Global.Phase.Oil:
                    if (thisBlock && upstream == block && diff_variable != Global.Variable.Pressure
                        ||
                        !thisBlock && upstream == neighbor && diff_variable != Global.Variable.Pressure)
                    {
                        if (diff_variable == Global.Variable.SaturationGas)
                        {
                            Kr = data.scal.GetKro(upstream.Sg[2], upstream.Sw[1], data.pvt.connateWaterSaturation);
                        }
                        else // saturation water
                        {
                            Kr = data.scal.GetKro(upstream.Sg[1], upstream.Sw[2], data.pvt.connateWaterSaturation);
                        }
                    }
                    else
                    {
                        Kr = upstream.Kro[1];
                    }
                    break;
                case Global.Phase.Gas:
                    if (thisBlock && upstream == block && diff_variable == Global.Variable.SaturationGas
                        ||
                        !thisBlock && upstream == neighbor && diff_variable == Global.Variable.SaturationGas)
                    {
                        Kr = upstream.Krg[2];
                    }
                    else
                    {
                        Kr = upstream.Krg[1];
                    }
                    break;
                default:
                    Kr = 1;
                    break;
            }

            return Kr;
        }

        public static double GetOilPotentialDifference(BaseBlock neighbor, BaseBlock block, double delta_z, double gravity)
        {
            return neighbor.P[1] - block.P[1] - delta_z * gravity;
        }

        public static double GetGasPotentialDifference(BaseBlock neighbor, BaseBlock block, double delta_z, double gravity)
        {
            return neighbor.GetPg(1, 1) - block.GetPg(1, 1) - delta_z * gravity;
        }

        public static double GetWaterPotentialDifference(BaseBlock neighbor, BaseBlock block, double delta_z, double gravity)
        {
            return neighbor.GetPw(1, 1) - block.GetPw(1, 1) - delta_z * gravity;
        }
    }
}
