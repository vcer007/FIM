using FIM.Core;

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
            // Note : Kro, Krw, Krg are all calculated based on Sg.

            BaseBlock upstream_block, downstream_block, neighbor_block;
            double transmissibility;

            double Kr, B, viscosity, Rso;
            double R = 0;

            double temp, accumulation_term;


            if (phase == Global.Phase.Oil)
            {
                for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
                {
                    if (block.neighbourBlocksIndices[i] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighbourBlocksIndices[i]];
                    transmissibility = block.transmissibility_list[i];

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

                    Kr = upstream_block.Kro[1];
                    B = 0.5 * (block.Bo[1] + neighbor_block.Bo[1]);
                    viscosity = 0.5 * (block.viscosityOil[1] + neighbor_block.viscosityOil[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1]);
                    block.transmissibility_terms_oil[i] = temp;

                    R += temp;

                }

                accumulation_term = 1 / (Global.a * data.timeStep) * ((block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                R -= accumulation_term;
            }
            else if (phase == Global.Phase.Water)
            {
                for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
                {
                    if (block.neighbourBlocksIndices[i] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighbourBlocksIndices[i]];
                    transmissibility = block.transmissibility_list[i];

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

                    Kr = upstream_block.Krw[1];
                    B = 0.5 * (block.Bw[1] + neighbor_block.Bw[1]);
                    viscosity = 0.5 * (block.viscosityWater[1] + neighbor_block.viscosityWater[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1]);
                    block.transmissibility_terms_water[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / (Global.a * data.timeStep) * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                R -= accumulation_term;
            }
            else
            {
                for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
                {
                    if (block.neighbourBlocksIndices[i] < 0)
                    {
                        continue;
                    }
                    neighbor_block = data.grid[block.neighbourBlocksIndices[i]];
                    transmissibility = block.transmissibility_list[i];

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

                    Kr = upstream_block.Krg[1];
                    B = 0.5 * (block.Bg[1] + neighbor_block.Bg[1]);
                    viscosity = 0.5 * (block.viscosityGas[1] + neighbor_block.viscosityGas[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1]);

                    if (data.solubleGasPresent)
                    {
                        Kr = upstream_block.Kro[1];
                        B = 0.5 * (block.Bo[1] + neighbor_block.Bo[1]);
                        viscosity = 0.5 * (block.viscosityOil[1] + neighbor_block.viscosityOil[1]);
                        Rso = 0.5 * (block.Rso[1] + neighbor_block.Rso[1]);

                        temp += Rso * transmissibility * Kr / (viscosity * B) * (neighbor_block.P[1] - block.P[1]);
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
            double R_plus = 0;

            double transmissibility_temp = 0;

            double transmissibility_term = 0;
            double accumulation_term = 0;

            double kr, viscosity, B, Rso;
            BaseBlock upstream_block, downstream_block, neighbor_block;

            double transmissiblity;

            #region Variable Dependencies
            int pd = 1, npd = 1, sod = 1, nsod = 1, swd = 1, nswd = 1, sgd = 1, nsgd = 1, swgd = 1, nswgd = 1;
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
                    swgd = 2;
                }
                else
                {
                    nsgd = 2;
                    nswgd = 2;
                }
            }
            else if (variable == Global.Variable.SaturationWater)
            {
                if (this_block)
                {
                    swd = 2;
                    swgd = 2;
                }
                else
                {
                    nswd = 2;
                    nswgd = 2;
                }
            }
            #endregion


            if (equation_phase == Global.Phase.Oil)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
                    {
                        if (block.neighbourBlocksIndices[i] == -1)
                        {
                            continue;
                        }

                        neighbor_block = data.grid[block.neighbourBlocksIndices[i]];

                        if (block.P[1] >= neighbor_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbor_block;

                            kr = upstream_block.Kro[sgd];
                        }
                        else
                        {
                            upstream_block = neighbor_block;
                            downstream_block = block;

                            kr = upstream_block.Kro[nsgd];
                        }

                        B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.timeStep * Global.a) * ((block.Vp[pd] * (1 - block.Sg[sgd] - block.Sw[swd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                else
                {
                    // transmissibility term
                    neighbor_block = data.grid[block.neighbourBlocksIndices[neighbour_block_index]];

                    if (block.P[1] >= neighbor_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbor_block;

                        kr = upstream_block.Kro[sgd];
                    }
                    else
                    {
                        upstream_block = neighbor_block;
                        downstream_block = block;

                        kr = upstream_block.Kro[nsgd];
                    }

                    B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                    viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
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
                    for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
                    {
                        if (block.neighbourBlocksIndices[i] == -1)
                        {
                            continue;
                        }

                        neighbor_block = data.grid[block.neighbourBlocksIndices[i]];

                        if (block.P[1] >= neighbor_block.P[1])
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

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);

                        if (data.solubleGasPresent)
                        {
                            if (upstream_block.index == block.index)
                            {
                                kr = upstream_block.Kro[sgd];
                            }
                            else
                            {
                                kr = upstream_block.Kro[nsgd];
                            }

                            Rso = 0.5 * (block.Rso[pd] + neighbor_block.Rso[npd]);
                            B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                            viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);

                            transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);
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
                    neighbor_block = data.grid[block.neighbourBlocksIndices[neighbour_block_index]];

                    if (block.P[1] >= neighbor_block.P[1])
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

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);

                    if (data.solubleGasPresent)
                    {
                        if (upstream_block.index == block.index)
                        {
                            kr = upstream_block.Kro[sgd];
                        }
                        else
                        {
                            kr = upstream_block.Kro[nsgd];
                        }

                        Rso = 0.5 * (block.Rso[pd] + neighbor_block.Rso[npd]);
                        B = 0.5 * (block.Bo[pd] + neighbor_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosityOil[pd] + neighbor_block.viscosityOil[npd]);

                        transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);
                    }

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
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
                    for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
                    {
                        if (block.neighbourBlocksIndices[i] == -1)
                        {
                            continue;
                        }

                        neighbor_block = data.grid[block.neighbourBlocksIndices[i]];

                        if (block.P[1] >= neighbor_block.P[1])
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

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.timeStep * Global.a) * ((block.Vp[pd] * block.Sw[swd] / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));

                }
                else
                {
                    // transmissibility term
                    neighbor_block = data.grid[block.neighbourBlocksIndices[neighbour_block_index]];

                    if (block.P[1] >= neighbor_block.P[1])
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

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbor_block.P[npd] - block.P[pd]);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighbourBlocksIndices.Length; i++)
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

                minusR[counter] = -1 * block.CalculateR(data, Global.Phase.Oil);
                minusR[counter + 1] = -1 * block.CalculateR(data, Global.Phase.Gas);
                minusR[counter + 2] = -1 * block.CalculateR(data, Global.Phase.Water);

                counter += data.phases.Length;
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
        public static void CalculateJacobi_Matrix(SimulationData data, double[] minusR, double[][] Jacobi)
        {
            // this is necessary as the direct solver used modifies the jacobi matrix.
            for (int i = 0; i < Jacobi.Length; i++)
            {
                for (int j = 0; j < Jacobi[i].Length; j++)
                {
                    Jacobi[i][j] = 0;
                }
            }

            BaseBlock block;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                #region Oil
                // with respect to P
                Jacobi[counter][data.phases.Length * block.index] = (block.Perturb(data, Global.Phase.Oil, -1, Global.Variable.Pressure) + minusR[counter]) / Global.EPSILON;
                // with respect to Sg
                Jacobi[counter][data.phases.Length * block.index + 1] = (block.Perturb(data, Global.Phase.Oil, -1, Global.Variable.SaturationGas) + minusR[counter]) / Global.EPSILON;
                // with respect to Sw
                Jacobi[counter][data.phases.Length * block.index + 2] = (block.Perturb(data, Global.Phase.Oil, -1, Global.Variable.SaturationWater) + minusR[counter]) / Global.EPSILON;

                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    if (block.neighbourBlocksIndices[j] >= 0)
                    {
                        // with respect to P
                        Jacobi[counter][data.phases.Length * block.neighbourBlocksIndices[j]] = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.Pressure) + minusR[counter]) / Global.EPSILON;
                        // with respect to Sg
                        Jacobi[counter][data.phases.Length * block.neighbourBlocksIndices[j] + 1] = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.SaturationGas) + minusR[counter]) / Global.EPSILON;
                        // with respect to Sw
                        Jacobi[counter][data.phases.Length * block.neighbourBlocksIndices[j] + 2] = (block.Perturb(data, Global.Phase.Oil, j, Global.Variable.SaturationWater) + minusR[counter]) / Global.EPSILON;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                Jacobi[counter + 1][data.phases.Length * block.index] = (block.Perturb(data, Global.Phase.Gas, -1, Global.Variable.Pressure) + minusR[counter + 1]) / Global.EPSILON;
                // with respect to Sg
                Jacobi[counter + 1][data.phases.Length * block.index + 1] = (block.Perturb(data, Global.Phase.Gas, -1, Global.Variable.SaturationGas) + minusR[counter + 1]) / Global.EPSILON;
                // with respect to Sw
                Jacobi[counter + 1][data.phases.Length * block.index + 2] = (block.Perturb(data, Global.Phase.Gas, -1, Global.Variable.SaturationWater) + minusR[counter + 1]) / Global.EPSILON;

                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    if (block.neighbourBlocksIndices[j] >= 0)
                    {
                        // with respect to P
                        Jacobi[counter + 1][data.phases.Length * block.neighbourBlocksIndices[j]] = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.Pressure) + minusR[counter + 1]) / Global.EPSILON;
                        // with respect to Sg
                        Jacobi[counter + 1][data.phases.Length * block.neighbourBlocksIndices[j] + 1] = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.SaturationGas) + minusR[counter + 1]) / Global.EPSILON;
                        // with respect to Sw
                        Jacobi[counter + 1][data.phases.Length * block.neighbourBlocksIndices[j] + 2] = (block.Perturb(data, Global.Phase.Gas, j, Global.Variable.SaturationWater) + minusR[counter + 1]) / Global.EPSILON;
                    }
                }
                #endregion
                #region Water
                // with respect to P
                Jacobi[counter + 2][data.phases.Length * block.index] = (block.Perturb(data, Global.Phase.Water, -1, Global.Variable.Pressure) + minusR[counter + 2]) / Global.EPSILON;
                // with respect to Sg
                Jacobi[counter + 2][data.phases.Length * block.index + 1] = (block.Perturb(data, Global.Phase.Water, -1, Global.Variable.SaturationGas) + minusR[counter + 2]) / Global.EPSILON;
                // with respect to Sw
                Jacobi[counter + 2][data.phases.Length * block.index + 2] = (block.Perturb(data, Global.Phase.Water, -1, Global.Variable.SaturationWater) + minusR[counter + 2]) / Global.EPSILON;

                for (int j = 0; j < block.neighbourBlocksIndices.Length; j++)
                {
                    if (block.neighbourBlocksIndices[j] >= 0)
                    {
                        // with respect to P
                        Jacobi[counter + 2][data.phases.Length * block.neighbourBlocksIndices[j]] = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.Pressure) + minusR[counter + 2]) / Global.EPSILON;
                        // with respect to Sg
                        Jacobi[counter + 2][data.phases.Length * block.neighbourBlocksIndices[j] + 1] = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.SaturationGas) + minusR[counter + 2]) / Global.EPSILON;
                        // with respect to Sw
                        Jacobi[counter + 2][data.phases.Length * block.neighbourBlocksIndices[j] + 2] = (block.Perturb(data, Global.Phase.Water, j, Global.Variable.SaturationWater) + minusR[counter + 2]) / Global.EPSILON;
                    }
                }
                #endregion

                counter += data.phases.Length;
            }

            //To-Do : add well production/injection derivatives.

        }

    }
}
