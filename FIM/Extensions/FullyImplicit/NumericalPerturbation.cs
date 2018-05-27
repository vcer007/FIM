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

            double Kr, B, viscosity, Rso;
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

                    if (GravityCalculation)
                    {
                        gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Oil, block.P[0]) + data.pvt.GetDensity(Global.Phase.Oil, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

                    transmissibility = block.transmissibility_list[i];

                    if (neighbor_block.P[1] - block.P[1] - delta_z * gravity <= 0)
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

                    R += temp;

                }

                accumulation_term = 1 / (Global.a * data.timeStep) * ((block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                R -= accumulation_term;
            }
            else if (phase == Global.Phase.Water)
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
                        gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Water, block.P[0]) + data.pvt.GetDensity(Global.Phase.Water, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

                    transmissibility = block.transmissibility_list[i];

                    if (neighbor_block.P[1] - block.P[1] - delta_z * gravity <= 0)
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
            else
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
                        gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Gas, block.P[0]) + data.pvt.GetDensity(Global.Phase.Gas, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

                    transmissibility = block.transmissibility_list[i];

                    if (neighbor_block.GetPg(1, 1) - block.GetPg(1, 1) - delta_z * gravity <= 0)
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
                        Kr = upstream_block.Kro[1];
                        B = 0.5 * (block.Bo[1] + neighbor_block.Bo[1]);
                        viscosity = 0.5 * (block.viscosityOil[1] + neighbor_block.viscosityOil[1]);
                        Rso = 0.5 * (block.Rso[1] + neighbor_block.Rso[1]);

                        if (GravityCalculation)
                        {
                            gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Oil, block.P[0]) + data.pvt.GetDensity(Global.Phase.Oil, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                        }

                        temp += Rso * transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1] - delta_z * gravity);
                        //temp += Rso * block.transmissibility_terms_oil[i];
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

            double kr = 1, viscosity, B, Rso;
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

                        if (GravityCalculation)
                        {
                            gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Oil, block.P[0]) + data.pvt.GetDensity(Global.Phase.Oil, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                        }
                        double delta_z = neighbor_block.Depth - block.Depth;

                        if (neighbor_block.P[1] - block.P[1] - delta_z * gravity <= 0)
                        {
                            upstream_block = block;
                            downstream_block = neighbor_block;

                            //kr = upstream_block.Kro[sgd];
                            if (variable == Global.Variable.SaturationGas)
                            {
                                if (sgd == 1)
                                {
                                    kr = upstream_block.Kro[sgd];
                                }
                                else if (sgd == 2)
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                                }
                            }
                            else if (variable == Global.Variable.SaturationWater)
                            {
                                if (swd == 1)
                                {
                                    kr = upstream_block.Kro[swd];
                                }
                                else if (swd == 2)
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                                }
                            }
                            else
                            {
                                kr = upstream_block.Kro[1];
                            }
                            
                        }
                        else
                        {
                            upstream_block = neighbor_block;
                            downstream_block = block;

                            //kr = upstream_block.Kro[nsgd];
                            if (variable == Global.Variable.SaturationGas)
                            {
                                if (nsgd == 1)
                                {
                                    kr = upstream_block.Kro[nsgd];
                                }
                                else
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                                }
                            }
                            else if (variable == Global.Variable.SaturationWater)
                            {
                                if (nsgd == 1)
                                {
                                    kr = upstream_block.Kro[nsgd];
                                }
                                else
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                                }
                            }
                            else
                            {
                                kr = upstream_block.Kro[1];
                            }
                        }

                        B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd] - delta_z * gravity);

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.timeStep * Global.a) * ((block.Vp[pd] * (1 - block.Sg[sgd] - block.Sw[swd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                else
                {
                    // transmissibility term
                    neighbor_block = data.grid[block.neighborBlocksIndices[neighbour_block_index]];

                    if (GravityCalculation)
                    {
                        gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Oil, block.P[0]) + data.pvt.GetDensity(Global.Phase.Oil, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (neighbor_block.P[1] - block.P[1] - delta_z * gravity <= 0)
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;

                        //kr = upstream_block.Kro[sgd];
                        if (variable == Global.Variable.SaturationGas)
                        {
                            if (sgd == 1)
                            {
                                kr = upstream_block.Kro[sgd];
                            }
                            else if (sgd == 2)
                            {
                                kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                            }
                        }
                        else if (variable == Global.Variable.SaturationWater)
                        {
                            if (swd == 1)
                            {
                                kr = upstream_block.Kro[swd];
                            }
                            else if (swd == 2)
                            {
                                kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                            }
                        }
                        else
                        {
                            kr = upstream_block.Kro[1];
                        }
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;

                        //kr = upstream_block.Kro[nsgd];
                        if (variable == Global.Variable.SaturationGas)
                        {
                            if (nsgd == 1)
                            {
                                kr = upstream_block.Kro[nsgd];
                            }
                            else
                            {
                                kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                            }
                        }
                        else if (variable == Global.Variable.SaturationWater)
                        {
                            if (nsgd == 1)
                            {
                                kr = upstream_block.Kro[nsgd];
                            }
                            else
                            {
                                kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                            }
                        }
                        else
                        {
                            kr = upstream_block.Kro[1];
                        }
                    }

                    B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                    viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd] - delta_z * gravity);

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

                        if (GravityCalculation)
                        {
                            gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Gas, block.P[0]) + data.pvt.GetDensity(Global.Phase.Gas, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                        }
                        double delta_z = neighbor_block.Depth - block.Depth;

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

                        B = 0.5 * (block.Bg[pd] + neighbor_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosityGas[pd] + neighbor_block.viscosityGas[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPg(npd, nsgd) - block.GetPg(pd, sgd) - delta_z * gravity);

                        if (data.solubleGasPresent)
                        {
                            if (upstream_block.index == block.index)
                            {
                                //kr = upstream_block.Kro[sgd];
                                if (variable == Global.Variable.SaturationGas)
                                {
                                    if (sgd == 1)
                                    {
                                        kr = upstream_block.Kro[sgd];
                                    }
                                    else if (sgd == 2)
                                    {
                                        kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                                    }
                                }
                                else if (variable == Global.Variable.SaturationWater)
                                {
                                    if (swd == 1)
                                    {
                                        kr = upstream_block.Kro[swd];
                                    }
                                    else if (swd == 2)
                                    {
                                        kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                                    }
                                }
                                else
                                {
                                    kr = upstream_block.Kro[1];
                                }
                            }
                            else
                            {
                                //kr = upstream_block.Kro[nsgd];
                                if (variable == Global.Variable.SaturationGas)
                                {
                                    if (nsgd == 1)
                                    {
                                        kr = upstream_block.Kro[nsgd];
                                    }
                                    else
                                    {
                                        kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                                    }
                                }
                                else if (variable == Global.Variable.SaturationWater)
                                {
                                    if (nsgd == 1)
                                    {
                                        kr = upstream_block.Kro[nsgd];
                                    }
                                    else
                                    {
                                        kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                                    }
                                }
                                else
                                {
                                    kr = upstream_block.Kro[1];
                                }
                            }

                            Rso = 0.5 * (block.Rso[pd] + neighbor_block.Rso[npd]);
                            B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                            viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);

                            if (GravityCalculation)
                            {
                                gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Oil, block.P[0]) + data.pvt.GetDensity(Global.Phase.Oil, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                            }

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

                    if (GravityCalculation)
                    {
                        gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Gas, block.P[0]) + data.pvt.GetDensity(Global.Phase.Gas, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

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

                    B = 0.5 * (block.Bg[pd] + neighbor_block.Bg[npd]);
                    viscosity = 0.5 * (block.viscosityGas[pd] + neighbor_block.viscosityGas[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.GetPg(npd, nsgd) - block.GetPg(pd, sgd) - delta_z * gravity);

                    if (data.solubleGasPresent)
                    {
                        if (upstream_block.index == block.index)
                        {
                            //kr = upstream_block.Kro[sgd];
                            if (variable == Global.Variable.SaturationGas)
                            {
                                if (sgd == 1)
                                {
                                    kr = upstream_block.Kro[sgd];
                                }
                                else if (sgd == 2)
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                                }
                            }
                            else if (variable == Global.Variable.SaturationWater)
                            {
                                if (swd == 1)
                                {
                                    kr = upstream_block.Kro[swd];
                                }
                                else if (swd == 2)
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                                }
                            }
                            else
                            {
                                kr = upstream_block.Kro[1];
                            }
                        }
                        else
                        {
                            //kr = upstream_block.Kro[nsgd];
                            if (variable == Global.Variable.SaturationGas)
                            {
                                if (nsgd == 1)
                                {
                                    kr = upstream_block.Kro[nsgd];
                                }
                                else
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[2], upstream_block.Sw[1], 0.12);
                                }
                            }
                            else if (variable == Global.Variable.SaturationWater)
                            {
                                if (nsgd == 1)
                                {
                                    kr = upstream_block.Kro[nsgd];
                                }
                                else
                                {
                                    kr = data.scal.GetKro(upstream_block.Sg[1], upstream_block.Sw[2], 0.12);
                                }
                            }
                            else
                            {
                                kr = upstream_block.Kro[1];
                            }
                        }

                        Rso = 0.5 * (block.Rso[pd] + neighbor_block.Rso[npd]);
                        B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);

                        if (GravityCalculation)
                        {
                            gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Oil, block.P[0]) + data.pvt.GetDensity(Global.Phase.Oil, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                        }

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

                        if (GravityCalculation)
                        {
                            gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Water, block.P[0]) + data.pvt.GetDensity(Global.Phase.Water, neighbor_block.P[0])) * Global.gamma_c * Global.alpha;
                        }
                        double delta_z = neighbor_block.Depth - block.Depth;

                        if (neighbor_block.GetPw(1, 1) - block.GetPw(1, 1) - delta_z * gravity <= 0)
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

                    if (GravityCalculation)
                    {
                        gravity = 0.5 * (data.pvt.GetDensity(Global.Phase.Water, block.GetPw(0, 0)) + data.pvt.GetDensity(Global.Phase.Water, neighbor_block.GetPw(0, 0))) * Global.gamma_c * Global.alpha;
                    }
                    double delta_z = neighbor_block.Depth - block.Depth;

                    if (neighbor_block.GetPw(1, 1) - block.GetPw(1, 1) - delta_z * gravity <= 0)
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

            //Jacobi.Clear();

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
                        temp = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.SaturationGas) + minusR[counter]) / Global.EPSILON_S;
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
                        temp = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.SaturationGas) + minusR[counter + 1]) / Global.EPSILON_S;
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
                        temp = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.SaturationGas) + minusR[counter + 2]) / Global.EPSILON_S;
                        Jacobi[counter + 2, place + 1] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;

                        // with respect to Sw
                        temp = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.SaturationWater) + minusR[counter + 2]) / Global.EPSILON_S;
                        Jacobi[counter + 2, place + 2] = temp > -Global.MINIMUM && temp < Global.MINIMUM ? Global.MINIMUM : temp;
                    }
                }
                #endregion
            }
        }

    }
}
