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

        public double[] area_list, delta_x_list, boundary_length_list;

        public double h, bulk_volume;

        public double[] Vp;

        // Well data
        public Global.BlockType type;
        public Global.WellType well_type;

        public double r_equivalent, WI;

        public double specified_flow_rate, specified_BHP;

        public double[] BHP;

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
            this.type = Global.BlockType.Normal_Block;
            this.well_type = Global.WellType.ShutIn;
            this.BHP = new double[steps_memory];
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
            // for the Sw value is constant.
            this.Sw[0] = Sw; this.Sw[1] = Sw; this.Sw[2] = Sw;
            this.Sg[0] = Sg; this.Sg[1] = Sg; this.Sg[2] = Sg + Global.epsilon;

            // Kro is only dependent on Sg. Sw is irreducible and won't increase.
            // if So is at max "Sg = 0" , then So can not increment up.
            this.Kro[0] = kr.getKr(Global.Phase.Oil, Sg); this.Kro[1] = this.Kro[0]; this.Kro[2] = kr.getKr(Global.Phase.Oil, this.Sg[2]);
            // Krw value is always equal to zero in this case.
            this.Krw[0] = 0; this.Krw[1] = this.Krw[0]; this.Krw[2] = this.Krw[0];
            this.Krg[0] = kr.getKr(Global.Phase.Gas, Sg); this.Krg[1] = this.Krg[0]; this.Krg[2] = kr.getKr(Global.Phase.Gas, this.Sg[2]);

            //this.Po = new double[steps_memory]; this.Pg = new double[steps_memory]; this.Pw = new double[steps_memory];

            // volumetric
            this.Vp[0] = this.bulk_volume * this.porosity[0]; this.Vp[1] = this.Vp[0]; this.Vp[2] = this.bulk_volume * this.porosity[2];

            // well
            this.BHP[0] = Well.WellData.calculatePwf(this, this.P[0], this.Kro[0], this.viscosity_oil[0]);
            this.BHP[1] = this.BHP[0];
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
            this.BHP[1] = Well.WellData.calculatePwf(this, this.P[1], this.Kro[1], this.viscosity_oil[1]);
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

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.So[1] / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                block.production_term_oil = block.q_oil[1];

                R -= accumulation_term;
                R -= block.production_term_oil;
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

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                block.q_water[1] = 0;
                block.production_term_water = block.q_water[1];


                R -= accumulation_term;
                R -= block.production_term_water;
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

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));

                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    accumulation_term += 1 / (Global.a * time_step) * ((block.Rso[1] * block.Vp[1] * block.So[1] / block.Bo[1]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                }

                if (block.type == Global.BlockType.Well_Block)
                {
                    if (block.well_type == Global.WellType.Production)
                    {
                        block.BHP[1] = Well.WellData.calculatePwf(block, block.P[1], block.Kro[1], block.viscosity_oil[1]);
                        block.q_gas[1] = Well.WellData.calculateFlow_Rate(block.P[1], block.BHP[1], block.Krg[1], block.viscosity_gas[1], block.WI);

                        // check for presence of soluble_gas in simulation_data
                        if (solubleGasPresent)
                        {
                            block.q_gas[1] += block.Rso[1] * block.q_oil[1];
                        }
                    }
                    else if (block.well_type == Global.WellType.Injection)
                    {

                    }
                }
                

                block.accumulation_term_gas = accumulation_term;
                block.production_term_gas = block.q_gas[1];

                R -= accumulation_term;
                R -= block.production_term_gas;
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


            #region pressure-dependent saturation-dependent, neighbour block pressure-dependent saturation-dependent.
            int pd = 1, npd = 1, sod = 1, nsod = 1, sgd = 1, nsgd = 1, swd = 1, nswd = 1;
            bool this_block = false;
            this_block = neighbour_block_index >= 0 ? false : true;

            if (variable == Global.Variable.Pressure)
            {
                if (!this_block)
                {
                    npd = 2;
                }
                else
                {
                    pd = 2;
                }
            }
            else if (variable == Global.Variable.Saturation_Oil)
            {
                if (!this_block)
                {
                    nsod = 2;
                }
                else
                {
                    sod = 2;
                }
            }
            else if(variable == Global.Variable.Saturation_Gas)
            {
                if (!this_block)
                {
                    nsgd = 2;
                }
                else
                {
                    sgd = 2;
                }
            }
            else
            {
                if (!this_block)
                {
                    nswd = 2;
                }
                else
                {
                    swd = 2;
                }
            }
            #endregion

            #region Oil
            if (phase == Global.Phase.Oil)
            {
                #region Block
                if (this_block)
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

                        Kr = upstream_block.Kro[sgd];
                        B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                        R += temp;

                    }

                    accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * (1 - block.Sw[1] - block.Sg[sgd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0])) + block.q_oil[1];

                    R -= accumulation_term;

                }
                #endregion
                #region NextBlock
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

                    Kr = upstream_block.Kro[sgd];
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
                #endregion
            }
            #endregion
            #region Water
            else if (phase == Global.Phase.Water)
            {
                #region Block
                if (this_block)
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

                        Kr = upstream_block.Krw[swd];
                        B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                        viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                        R += temp;
                    }

                    // here there is no Sw "instead (1 - So - Sg)". so when differentiating with respect to either So or Sg, we increment either one of them.
                    accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * (1 - block.So[sod] - block.Sg[1]) / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0])) + block.q_water[1];

                    R -= accumulation_term;
                }
                #endregion
                #region NextBlock
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

                    Kr = upstream_block.Krw[swd];
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
                #endregion
            }
            #endregion
            #region Gas
            else
            {
                #region Block
                if (this_block)
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

                        Kr = upstream_block.Krg[sgd];
                        B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (block.P[pd] - neighbour_block.P[npd]);

                        R += temp;
                    }

                    accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * block.Sg[sgd] / block.Bg[pd]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                    // check for presence of soluble_gas in simulation_data
                    if (solubleGasPresent)
                    {
                        accumulation_term += 1 / (Global.a * time_step) * ((block.Rso[pd] * block.Vp[pd] * (1 - block.Sw[1] - block.Sg[sgd]) / block.Bo[pd]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                    }

                    if (block.type == Global.BlockType.Well_Block)
                    {
                        if (block.well_type == Global.WellType.Production)
                        {
                            block.BHP[2] = Well.WellData.calculatePwf(block, block.P[pd], block.Kro[sgd], block.viscosity_oil[pd]);
                            block.q_gas[2] = Well.WellData.calculateFlow_Rate(block.P[pd], block.BHP[2], block.Krg[sgd], block.viscosity_gas[pd], block.WI);

                            // check for presence of soluble_gas in simulation_data
                            if (solubleGasPresent)
                            {
                                block.q_gas[2] += block.Rso[pd] * block.q_oil[1];
                            }
                        }
                        else if (block.well_type == Global.WellType.Injection)
                        {

                        }

                        production_term_gas = block.q_gas[2];
                    }
                    

                    R -= accumulation_term;
                    R -= production_term_gas;
                }
                #endregion
                #region NextBlock
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

                    Kr = upstream_block.Krg[sgd];
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
                #endregion
            }
            #endregion

            return R;

        }
    }
}
