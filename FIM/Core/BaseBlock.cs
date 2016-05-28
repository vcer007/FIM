using FIM.Fluid;
using FIM.Rock;

namespace FIM.Core
{
    /// <summary>
    /// The base class that contains all the common properties that any block should have "regardless of specific block type".
    /// </summary>
    class BaseBlock
    {
        /// <summary>
        /// A list of the indices of neighbouring blocks.
        /// </summary>
        public int index;
        public int layer;
        public int[] neighbour_blocks_indices;


        // Petro-physical properties.


        // Rock properties

        public double[] porosity;

        public double[] permeability_list;

        // Fluid properties

        public double[] Sw, So, Sg;

        public double[] Krw, Kro, Krg;

        public double[] Bw, Bo, Bg;

        public double[] viscosity_water, viscosity_oil, viscosity_gas;

        public double[] Rso;

        public double[] P, Pw, Po, Pg;


        // Volumetric data

        public double[] area_list, delta_x_list;

        public double h, bulk_volume;

        public double[] Vp;

        // Well data

        public double[] flow_rate, BHP, J;

        public double well_radius, skin;

        public double[] q_water, q_oil, q_gas;

        // list of the blocks transmissibilities, brginning with the top then bottom ones.

        public double[] transmissibility_list;


        // these variables are used to save repeated calculations unnecessarily.
        public double[] transmissibility_terms_oil, transmissibility_terms_water, transmissibility_terms_gas;
        public double accumulation_term_oil, accumulation_term_water, accumulation_term_gas;
        public double production_term_oil, production_term_water, production_term_gas;


        // public class constructor initializing the array sizes.
        public BaseBlock()
        {
            // we store only three time steps "current time step "n0" and next time step "n1" or the first iteration" and the third is used for jacobians generation.
            const int steps_memory = 3;

            // rock
            this.porosity = new double[steps_memory];

            // fluid
            this.Bo = new double[steps_memory]; this.Bg = new double[steps_memory]; this.Bw = new double[steps_memory];
            this.So = new double[steps_memory]; this.Sg = new double[steps_memory]; this.Sw = new double[steps_memory];
            this.Kro = new double[steps_memory]; this.Krg = new double[steps_memory]; this.Krw = new double[steps_memory];
            this.viscosity_oil = new double[steps_memory]; this.viscosity_gas = new double[steps_memory]; this.viscosity_water = new double[steps_memory];
            this.Rso = new double[steps_memory];

            this.P = new double[steps_memory];  this.Po = new double[steps_memory]; this.Pg = new double[steps_memory]; this.Pw = new double[steps_memory];

            // volumetric
            this.Vp = new double[steps_memory];

            // well
            this.flow_rate = new double[steps_memory]; this.BHP = new double[steps_memory]; this.J = new double[steps_memory];
            this.q_oil = new double[steps_memory]; this.q_gas = new double[steps_memory]; this.q_water = new double[steps_memory];
        }

        public void updateProperties(PVT pvt, Kr kr, Porosity porosity, double P, double Sw, double So, double Sg)
        {
            // pressure
            this.P[0] = P; this.P[1] = P; this.P[2] = P + Global.epsilon;

            // rock
            this.porosity[0] = porosity.getPorosity(P); this.porosity[1] = this.porosity[0]; this.porosity[2] = porosity.getPorosity(this.P[2]);

            // fluid
            this.Bo[0] = pvt.getFVF(Global.Phase.Oil, P); this.Bo[1] = this.Bo[0]; this.Bo[2] = pvt.getFVF(Global.Phase.Oil, this.P[2]);
            this.Bw[0] = pvt.getFVF(Global.Phase.Water, P); this.Bw[1] = this.Bw[0]; this.Bw[2] = pvt.getFVF(Global.Phase.Water, this.P[2]);
            this.Bg[0] = pvt.getFVF(Global.Phase.Gas, P); this.Bg[1] = this.Bg[0]; this.Bg[2] = pvt.getFVF(Global.Phase.Gas, this.P[2]);

            this.viscosity_oil[0] = pvt.getViscosity(Global.Phase.Oil, P); this.viscosity_oil[1] = this.viscosity_oil[0]; this.viscosity_oil[2] = pvt.getViscosity(Global.Phase.Oil, this.P[2]);
            this.viscosity_water[0] = pvt.getViscosity(Global.Phase.Water, P); this.viscosity_water[1] = this.viscosity_water[0]; this.viscosity_water[2] = pvt.getViscosity(Global.Phase.Water, this.P[2]);
            this.viscosity_gas[0] = pvt.getViscosity(Global.Phase.Gas, P); this.viscosity_gas[1] = this.viscosity_gas[0]; this.viscosity_gas[2] = pvt.getViscosity(Global.Phase.Gas, this.P[2]);

            this.Rso[0] = pvt.getRs(Global.Phase.Oil, P); this.Rso[1] = this.Rso[0]; this.Rso[2] = pvt.getRs(Global.Phase.Oil, this.P[2]);

            this.So[0] = So; this.So[1] = So; this.So[2] = So + Global.epsilon;
            this.Sw[0] = Sw; this.Sw[1] = Sw; this.Sw[2] = Sw + Global.epsilon;
            this.Sg[0] = Sg; this.Sg[1] = Sg; this.Sg[2] = Sg + Global.epsilon;

            this.Kro[0] = kr.getKr(Global.Phase.Oil, Sg); this.Kro[1] = this.Kro[0]; this.Kro[2] = kr.getKr(Global.Phase.Oil, this.Sg[2]);
            this.Krw[0] = kr.getKr(Global.Phase.Water, Sg); this.Krw[1] = this.Krw[0]; this.Krw[2] = kr.getKr(Global.Phase.Water, this.Sg[2]);
            this.Krg[0] = kr.getKr(Global.Phase.Gas, Sg); this.Krg[1] = this.Krg[0]; this.Krg[2] = kr.getKr(Global.Phase.Gas, this.Sg[2]);

            //this.Po = new double[steps_memory]; this.Pg = new double[steps_memory]; this.Pw = new double[steps_memory];

            // volumetric
            this.Vp[0] = this.bulk_volume * this.porosity[0]; this.Vp[1] = this.Vp[0]; this.Vp[2] = this.bulk_volume * this.porosity[2];

            // well
        }

        public void updateProperties_n1(PVT pvt, Kr kr, Porosity porosity, double P, double Sw, double So, double Sg)
        {
            // pressure
            this.P[1] = P;

            // rock
            this.porosity[1] = porosity.getPorosity(P);

            // fluid
            this.Bo[1] = pvt.getFVF(Global.Phase.Oil, P);
            this.Bw[1] = pvt.getFVF(Global.Phase.Water, P);
            this.Bg[1] = pvt.getFVF(Global.Phase.Gas, P);

            this.viscosity_oil[1] = pvt.getViscosity(Global.Phase.Oil, P);
            this.viscosity_water[1] = pvt.getViscosity(Global.Phase.Water, P);
            this.viscosity_gas[1] = pvt.getViscosity(Global.Phase.Gas, P);

            this.Rso[1] = pvt.getRs(Global.Phase.Oil, P);

            this.So[1] = So;
            this.Sw[1] = Sw;
            this.Sg[1] = Sg;

            this.Kro[1] = kr.getKr(Global.Phase.Oil, Sg);
            this.Krw[1] = kr.getKr(Global.Phase.Water, Sg);
            this.Krg[1] = kr.getKr(Global.Phase.Gas, Sg);

            //this.Po = new double[steps_memory]; this.Pg = new double[steps_memory]; this.Pw = new double[steps_memory];

            // volumetric
            this.Vp[1] = this.bulk_volume * this.porosity[1];

            // well
        }

        public double calculateR(BaseBlock[] grid, Global.Phase phase, double time_step, bool solubleGasPresent = true)
        {
            // Note : Kro, Krw, Krg are all calculated based on Sg.

            BaseBlock block = this;

            BaseBlock upstream_block, downstream_block, neighbour_block;
            double transmissibility;

            double Kr, B, viscosity;
            double R = 0;

            double temp, accumulation_term, production_term;


            if (phase == Global.Phase.Oil)
            {
                for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                {
                    if (neighbour_blocks_indices[i] < 0)
                    {
                        continue;
                    }
                    neighbour_block = grid[block.neighbour_blocks_indices[i]];
                    transmissibility = block.transmissibility_list[i];

                    if (block.P[0] >= neighbour_block.P[0])
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

                    temp = transmissibility * Kr / (viscosity * B) * (block.P[1] - neighbour_block.P[1]);
                    block.transmissibility_terms_oil[i] = temp;

                    R += temp;

                }

                accumulation_term = 1 / time_step * ((block.Vp[1] * block.So[1] / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0])) + block.q_oil[1];
                block.accumulation_term_oil = accumulation_term;

                R -= accumulation_term;
            }
            else if (phase == Global.Phase.Water)
            {
                for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                {
                    if (neighbour_blocks_indices[i] < 0)
                    {
                        continue;
                    }
                    neighbour_block = grid[block.neighbour_blocks_indices[i]];
                    transmissibility = block.transmissibility_list[i];

                    if (block.P[0] >= neighbour_block.P[0])
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

                    temp = transmissibility * Kr / (viscosity * B) * (block.P[1] - neighbour_block.P[1]);
                    block.transmissibility_terms_water[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / time_step * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0])) + block.q_water[1];
                block.accumulation_term_water = accumulation_term;

                R -= accumulation_term;
            }
            else
            {
                for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                {
                    if (neighbour_blocks_indices[i] < 0)
                    {
                        continue;
                    }
                    neighbour_block = grid[block.neighbour_blocks_indices[i]];
                    transmissibility = block.transmissibility_list[i];

                    if (block.P[0] >= neighbour_block.P[0])
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

                    temp = transmissibility * Kr / (viscosity * B) * (block.P[1] - neighbour_block.P[1]);
                    block.transmissibility_terms_gas[i] = temp;

                    R += temp;
                }

                accumulation_term = 1 / time_step * ((block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Vp[0] * block.Sg[0] / block.Bg[0])) + block.q_gas[1];

                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    accumulation_term += 1 / time_step * ((block.Rso[1] * block.Vp[1] * block.So[1] / block.Bo[1]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                }

                block.accumulation_term_gas = accumulation_term;

                R -= accumulation_term;
            }

            return R;

        }

        public double calculateDerivative(BaseBlock[] grid, Global.Phase phase, double time_step, int neighbour_block_index, Global.Variable variable, bool solubleGasPresent = true)
        {
            BaseBlock block = this;

            BaseBlock upstream_block, downstream_block, neighbour_block;
            double transmissibility;

            double Kr, B, viscosity;
            double R = 0;

            double temp, accumulation_term, production_term;

            // pressure-dependent saturation-dependent, neighbour block pressure-dependent saturation-dependent.
            int pd = 1, sd = 1, npd = 1, nsd = 1;

            if (variable == Global.Variable.Pressure)
            {
                if (neighbour_block_index >= 0)
                {
                    npd = 2;
                }
                else
                {
                    pd = 2;
                }
            }
            else
            {
                if (neighbour_block_index >= 0)
                {
                    nsd = 2;
                }
                else
                {
                    sd = 2;
                }
            }
            

            #region Oil
            if (phase == Global.Phase.Oil)
            {
                if (pd == 2 || sd == 2)
                {
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (neighbour_blocks_indices[i] < 0)
                        {
                            continue;
                        }
                        neighbour_block = grid[block.neighbour_blocks_indices[i]];
                        transmissibility = block.transmissibility_list[i];

                        if (block.P[0] >= neighbour_block.P[0])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;
                        }

                        Kr = upstream_block.Kro[sd];
                        B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                        R += temp;

                    }

                    accumulation_term = 1 / time_step * ((block.Vp[pd] * block.So[sd] / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0])) + block.q_oil[1];

                    R -= accumulation_term;
                }
                else
                {
                    neighbour_block = grid[block.neighbour_blocks_indices[neighbour_block_index]];
                    transmissibility = block.transmissibility_list[neighbour_block_index];

                    if (block.P[0] >= neighbour_block.P[0])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Kro[sd];
                    B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                    viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                    temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                    R += temp;

                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i != neighbour_block_index)
                        {
                            R += block.transmissibility_terms_oil[i];
                        }
                    }

                    R -= block.accumulation_term_oil;
                    R -= block.production_term_oil;
                }
                
            }
            #endregion
            #region Water
            else if (phase == Global.Phase.Water)
            {
                if (pd == 2 || sd == 2)
                {
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (neighbour_blocks_indices[i] < 0)
                        {
                            continue;
                        }
                        neighbour_block = grid[block.neighbour_blocks_indices[i]];
                        transmissibility = block.transmissibility_list[i];

                        if (block.P[0] >= neighbour_block.P[0])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;
                        }

                        Kr = upstream_block.Krw[sd];
                        B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                        viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                        R += temp;
                    }

                    accumulation_term = 1 / time_step * ((block.Vp[pd] * block.Sw[sd] / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0])) + block.q_water[1];

                    R -= accumulation_term;
                }
                else
                {
                    neighbour_block = grid[block.neighbour_blocks_indices[neighbour_block_index]];
                    transmissibility = block.transmissibility_list[neighbour_block_index];

                    if (block.P[0] >= neighbour_block.P[0])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Krw[sd];
                    B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                    viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);

                    temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                    R += temp;

                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i != neighbour_block_index)
                        {
                            R += block.transmissibility_terms_water[i];
                        }
                    }

                    R -= block.accumulation_term_water;
                    R -= block.production_term_water;
                }
                
            }
            #endregion
            #region Gas
            else
            {
                if (pd == 2 || sd == 2)
                {
                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (neighbour_blocks_indices[i] < 0)
                        {
                            continue;
                        }
                        neighbour_block = grid[block.neighbour_blocks_indices[i]];
                        transmissibility = block.transmissibility_list[i];

                        if (block.P[0] >= neighbour_block.P[0])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;
                        }

                        Kr = upstream_block.Krg[sd];
                        B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                        R += temp;
                    }

                    accumulation_term = 1 / time_step * ((block.Vp[pd] * block.Sg[sd] / block.Bg[pd]) - (block.Vp[0] * block.Sg[0] / block.Bg[0])) + block.q_gas[1];

                    // check for presence of soluble_gas in simulation_data
                    if (solubleGasPresent)
                    {
                        accumulation_term += 1 / time_step * ((block.Rso[pd] * block.Vp[pd] * block.So[sd] / block.Bo[pd]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                    }

                    block.accumulation_term_gas = accumulation_term;

                    R -= accumulation_term;
                }
                else
                {
                    neighbour_block = grid[block.neighbour_blocks_indices[neighbour_block_index]];
                    transmissibility = block.transmissibility_list[neighbour_block_index];

                    if (block.P[0] >= neighbour_block.P[0])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;
                    }

                    Kr = upstream_block.Krg[sd];
                    B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                    viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);

                    temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                    R += temp;

                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i != neighbour_block_index)
                        {
                            R += block.transmissibility_terms_gas[i];
                        }
                    }

                    R -= block.accumulation_term_gas;
                    R -= block.production_term_gas;
                }
                
            }
            #endregion

            return R;

        }
    }
}
