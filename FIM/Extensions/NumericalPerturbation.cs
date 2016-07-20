using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace FIM.Extensions
{
    /// <summary>
    /// This class contains extension methods to the <see cref="BaseBlock"/> to add support for perturbation.
    /// </summary>
    public static class NumericalPerturbation
    {

        private static double calculateR(this BaseBlock block, SimulationData data, Global.Phase phase)
        {
            // Note : Kro, Krw, Krg are all calculated based on Sg.

            BaseBlock upstream_block, downstream_block, neighbour_block;

            double transmissibility;

            double Kr, B, viscosity, Rso;
            double R = 0;

            double temp, accumulation_term, production_term = 0;


            if (phase == Global.Phase.Oil)
            {
                for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                {
                    if (block.neighbour_blocks_indices[i] < 0)
                    {
                        continue;
                    }
                    neighbour_block = data.grid[block.neighbour_blocks_indices[i]];
                    transmissibility = block.transmissibility_list[i];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Kro[1];
                    B = 0.5 * (block.Bo[1] + neighbour_block.Bo[1]);
                    viscosity = 0.5 * (block.viscosity_oil[1] + neighbour_block.viscosity_oil[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[1] - block.P[1]);
                    block.transmissibility_terms_oil[i] = temp;

                    R += temp;

                }

                accumulation_term = 1 / (Global.a * data.time_step) * ((block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                //production_term = block.q_oil[1];
                block.production_term_oil = production_term;

                R -= (accumulation_term /*+ production_term*/);
            }
            else if (phase == Global.Phase.Water)
            {
                for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                {
                    if (block.neighbour_blocks_indices[i] < 0)
                    {
                        continue;
                    }
                    neighbour_block = data.grid[block.neighbour_blocks_indices[i]];
                    transmissibility = block.transmissibility_list[i];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Krw[1];
                    B = 0.5 * (block.Bw[1] + neighbour_block.Bw[1]);
                    viscosity = 0.5 * (block.viscosity_water[1] + neighbour_block.viscosity_water[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[1] - block.P[1]);
                    block.transmissibility_terms_water[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / (Global.a * data.time_step) * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                //production_term = block.q_water[1];
                block.production_term_water = production_term;


                R -= (accumulation_term /*+ production_term*/);
            }
            else
            {
                for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                {
                    if (block.neighbour_blocks_indices[i] < 0)
                    {
                        continue;
                    }
                    neighbour_block = data.grid[block.neighbour_blocks_indices[i]];
                    transmissibility = block.transmissibility_list[i];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Krg[1];
                    B = 0.5 * (block.Bg[1] + neighbour_block.Bg[1]);
                    viscosity = 0.5 * (block.viscosity_gas[1] + neighbour_block.viscosity_gas[1]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[1] - block.P[1]);

                    if (data.solubleGasPresent)
                    {
                        //Kr = upstream_block.Kro[1];
                        //B = 0.5 * (block.Bo[1] + neighbour_block.Bo[1]);
                        //viscosity = 0.5 * (block.viscosity_oil[1] + neighbour_block.viscosity_oil[1]);
                        Rso = 0.5 * (block.Rso[1] + neighbour_block.Rso[1]);

                        //temp += Rso * transmissibility * Kr / (viscosity * B) * (neighbour_block.P[1] - block.P[1]);
                        temp += Rso * block.transmissibility_terms_oil[i];
                    }

                    block.transmissibility_terms_gas[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / (Global.a * data.time_step) * ((block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                // check for presence of soluble_gas in simulation_data
                if (data.solubleGasPresent)
                {
                    accumulation_term += 1 / (Global.a * data.time_step) * ((block.Rso[1] * block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                block.accumulation_term_gas = accumulation_term;


                //if (block.type == Global.BlockType.Well_Block)
                //{
                //    if (block.well_type == Global.WellType.Production)
                //    {
                //        production_term = block.q_gas[1];
                //        // check for presence of soluble_gas in simulation_data
                //        if (data.solubleGasPresent)
                //        {
                //            production_term += block.Rso[1] * block.q_oil[1];
                //        }
                //    }
                //    else if (block.well_type == Global.WellType.Injection)
                //    {
                //        production_term = block.specified_flow_rate;
                //    }
                //}
                block.production_term_gas = production_term;


                R -= (accumulation_term /*+ production_term*/);
            }

            return R;

        }

        private static double perturb(this BaseBlock block, SimulationData data, Global.Phase phase, int neighbour_block_index, Global.Variable variable)
        {
            double R_plus = 0;

            double transmissibility_temp = 0;

            double transmissibility_term = 0;
            double accumulation_term = 0;
            double production_term = 0;

            double kr, viscosity, B, Rso;
            BaseBlock upstream_block, downstream_block, neighbour_block;

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
            else if (variable == Global.Variable.Saturation_Gas)
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
            else if (variable == Global.Variable.Saturation_Water)
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


            if (phase == Global.Phase.Oil)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (block.neighbour_blocks_indices[i] == -1)
                        {
                            continue;
                        }

                        neighbour_block = data.grid[block.neighbour_blocks_indices[i]];

                        if (block.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;

                            kr = upstream_block.Kro[sgd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;

                            kr = upstream_block.Kro[nsgd];
                        }

                        B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                        //block.transmissibility_terms_oil[i] = transmissibility_temp;

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.time_step * Global.a) * ((block.Vp[pd] * (1 - block.Sg[sgd] - block.Sw[swd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0]));

                    //block.accumulation_term_oil = accumulation_term;
                    // production term
                    //production_term = block.q_oil[1];

                    //block.production_term_oil = production_term;
                }
                else
                {
                    // transmissibility term
                    neighbour_block = data.grid[block.neighbour_blocks_indices[neighbour_block_index]];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;

                        kr = upstream_block.Kro[sgd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;

                        kr = upstream_block.Kro[nsgd];
                    }

                    B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                    viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += block.transmissibility_terms_oil[i];
                    }
                    // accumulation term
                    accumulation_term = block.accumulation_term_oil;
                    // production term
                    //production_term = block.production_term_oil;

                }
            }
            else if (phase == Global.Phase.Gas)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (block.neighbour_blocks_indices[i] == -1)
                        {
                            continue;
                        }

                        neighbour_block = data.grid[block.neighbour_blocks_indices[i]];

                        if (block.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;

                            kr = upstream_block.Krg[sgd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;

                            kr = upstream_block.Krg[nsgd];
                        }

                        B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

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

                            Rso = 0.5 * (block.Rso[pd] + neighbour_block.Rso[npd]);
                            B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                            viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                            transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);
                        }

                        //block.transmissibility_terms_gas[i] = transmissibility_temp;

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.time_step * Global.a) * ((block.Vp[pd] * block.Sg[sgd] / block.Bg[pd]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));

                    if (data.solubleGasPresent)
                    {
                        accumulation_term += 1 / (data.time_step * Global.a) * (block.Rso[pd] * (block.Vp[pd] * (1 - block.Sg[sgd] - block.Sw[swd]) / block.Bo[pd]) - block.Rso[0] * (block.Vp[0] * block.So[0] / block.Bo[0]));
                    }

                    //block.accumulation_term_gas = accumulation_term;
                    // production term
                    //if (block.type == Global.BlockType.Well_Block)
                    //{
                    //    if (block.well_type == Global.WellType.Production)
                    //    {
                    //        production_term = block.q_gas[1];

                    //        if (data.solubleGasPresent)
                    //        {
                    //            production_term += block.Rso[pd] * block.q_oil[1];
                    //        }
                    //    }
                    //    else
                    //    {
                    //        production_term = block.specified_flow_rate;
                    //    }
                    //}


                    //block.production_term_gas = production_term;
                }
                else
                {
                    // transmissibility term
                    neighbour_block = data.grid[block.neighbour_blocks_indices[neighbour_block_index]];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;

                        kr = upstream_block.Krg[sgd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;

                        kr = upstream_block.Krg[nsgd];
                    }

                    B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                    viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

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

                        Rso = 0.5 * (block.Rso[pd] + neighbour_block.Rso[npd]);
                        B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                        transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);
                    }

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += block.transmissibility_terms_gas[i];
                    }
                    // accumulation term
                    accumulation_term = block.accumulation_term_gas;
                    // production term
                    production_term = block.production_term_gas;
                }
            }
            else if (phase == Global.Phase.Water)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (block.neighbour_blocks_indices[i] == -1)
                        {
                            continue;
                        }

                        neighbour_block = data.grid[block.neighbour_blocks_indices[i]];

                        if (block.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;

                            kr = upstream_block.Krw[swd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;

                            kr = upstream_block.Krw[nswd];
                        }

                        B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                        viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);
                        transmissiblity = block.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                        //block.transmissibility_terms_water[i] = transmissibility_temp;

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.time_step * Global.a) * ((block.Vp[pd] * block.Sw[swd] / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));

                    //block.accumulation_term_water = accumulation_term;
                    // production term
                    //production_term = block.q_water[1];

                    //block.production_term_water = production_term;
                }
                else
                {
                    // transmissibility term
                    neighbour_block = data.grid[block.neighbour_blocks_indices[neighbour_block_index]];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;

                        kr = upstream_block.Krw[swd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;

                        kr = upstream_block.Krw[nswd];
                    }

                    B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                    viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);
                    transmissiblity = block.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += block.transmissibility_terms_water[i];
                    }
                    // accumulation term
                    accumulation_term = block.accumulation_term_water;
                    // production term
                    production_term = block.production_term_water;
                }
            }


            R_plus = transmissibility_term - accumulation_term /*- production_term*/;
            return R_plus;
        }


        public static void calculateMinusR_Matrix(SimulationData data, double[] minus_R)
        {
            //for (int i = 0; i < minus_R.Length; i++)
            //{
            //    minus_R[i] = 0;
            //}

            BaseBlock block;
            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                minus_R[counter] = -1 * block.calculateR(data, Global.Phase.Oil);
                minus_R[counter + 1] = -1 * block.calculateR(data, Global.Phase.Gas);
                minus_R[counter + 2] = -1 * block.calculateR(data, Global.Phase.Water);

                counter += data.phases.Length;
            }

            //To-Do : add well production/injection.
        }

        public static void calculateJacobi_Matrix(SimulationData data, double[] minus_R, double[][] jacobians)
        {
            //for (int i = 0; i < jacobians.Length; i++)
            //{
            //    for (int j = 0; j < jacobians[i].Length; j++)
            //    {
            //        jacobians[i][j] = 0;
            //    }
            //}
            //int size = data.grid.Length * data.phases.Length;
            ///*double[][] */jacobians = new double[size][];

            BaseBlock block;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                //jacobians[counter] = new double[size];
                //jacobians[counter + 1] = new double[size];
                //jacobians[counter + 2] = new double[size];

                #region Oil
                // with respect to P
                jacobians[counter][data.phases.Length * block.index] = (block.perturb(data, Global.Phase.Oil, -1, Global.Variable.Pressure) + minus_R[counter]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter][data.phases.Length * block.index + 1] = (block.perturb(data, Global.Phase.Oil, -1, Global.Variable.Saturation_Gas) + minus_R[counter]) / Global.epsilon;
                // with respect to Sw
                jacobians[counter][data.phases.Length * block.index + 2] = (block.perturb(data, Global.Phase.Oil, -1, Global.Variable.Saturation_Water) + minus_R[counter]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j]] = (block.perturb(data, Global.Phase.Oil, j, Global.Variable.Pressure) + minus_R[counter]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = (block.perturb(data, Global.Phase.Oil, j, Global.Variable.Saturation_Gas) + minus_R[counter]) / Global.epsilon;
                        // with respect to Sw
                        jacobians[counter][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.perturb(data, Global.Phase.Oil, j, Global.Variable.Saturation_Water) + minus_R[counter]) / Global.epsilon;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][data.phases.Length * block.index] = (block.perturb(data, Global.Phase.Gas, -1, Global.Variable.Pressure) + minus_R[counter + 1]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter + 1][data.phases.Length * block.index + 1] = (block.perturb(data, Global.Phase.Gas, -1, Global.Variable.Saturation_Gas) + minus_R[counter + 1]) / Global.epsilon;
                // with respect to Sw
                jacobians[counter + 1][data.phases.Length * block.index + 2] = (block.perturb(data, Global.Phase.Gas, -1, Global.Variable.Saturation_Water) + minus_R[counter + 1]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j]] = (block.perturb(data, Global.Phase.Gas, j, Global.Variable.Pressure) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = (block.perturb(data, Global.Phase.Gas, j, Global.Variable.Saturation_Gas) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sw
                        jacobians[counter + 1][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.perturb(data, Global.Phase.Gas, j, Global.Variable.Saturation_Water) + minus_R[counter + 1]) / Global.epsilon;
                    }
                }
                #endregion
                #region Water
                // with respect to P
                jacobians[counter + 2][data.phases.Length * block.index] = (block.perturb(data, Global.Phase.Water, -1, Global.Variable.Pressure) + minus_R[counter + 2]) / Global.epsilon;
                // with respect to Sg
                jacobians[counter + 2][data.phases.Length * block.index + 1] = (block.perturb(data, Global.Phase.Water, -1, Global.Variable.Saturation_Gas) + minus_R[counter + 2]) / Global.epsilon;
                // with respect to Sw
                jacobians[counter + 2][data.phases.Length * block.index + 2] = (block.perturb(data, Global.Phase.Water, -1, Global.Variable.Saturation_Water) + minus_R[counter + 2]) / Global.epsilon;

                for (int j = 0; j < block.neighbour_blocks_indices.Length; j++)
                {
                    if (block.neighbour_blocks_indices[j] >= 0)
                    {
                        // with respect to P
                        jacobians[counter + 2][data.phases.Length * block.neighbour_blocks_indices[j]] = (block.perturb(data, Global.Phase.Water, j, Global.Variable.Pressure) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sg
                        jacobians[counter + 2][data.phases.Length * block.neighbour_blocks_indices[j] + 1] = (block.perturb(data, Global.Phase.Water, j, Global.Variable.Saturation_Gas) + minus_R[counter + 1]) / Global.epsilon;
                        // with respect to Sw
                        jacobians[counter + 2][data.phases.Length * block.neighbour_blocks_indices[j] + 2] = (block.perturb(data, Global.Phase.Water, j, Global.Variable.Saturation_Water) + minus_R[counter + 1]) / Global.epsilon;
                    }
                }
                #endregion

                counter += data.phases.Length;
            }

            //return jacobians;
            //To-Do : add well production/injection derivatives.

        }

    }
}
