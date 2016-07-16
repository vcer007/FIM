using System;
using FIM.Fluid;
using FIM.Rock;

namespace FIM.Core
{
    /// <summary>
    /// The base class that contains all the common properties that any block should have "regardless of specific block type".
    /// </summary>
    public class BaseBlock
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
        public double P_previous_step;


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

        public double total_qg, GOR;

        // list of the blocks transmissibilities, brginning with the top then bottom ones.

        public double[] transmissibility_list;

        public bool saturated = false;


        // these variables are used to save repeated calculations unnecessarily.
        public double[] transmissibility_terms_oil, transmissibility_terms_water, transmissibility_terms_gas;
        public double accumulation_term_oil, accumulation_term_water, accumulation_term_gas;
        public double production_term_oil, production_term_water, production_term_gas;

        public double dp, dso, dsw, dsg;


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
            this.Rso = new double[steps_memory + 1];

            this.P = new double[steps_memory + 1];  this.Po = new double[steps_memory]; this.Pg = new double[steps_memory]; this.Pw = new double[steps_memory];

            // volumetric
            this.Vp = new double[steps_memory];

            // well
            this.type = Global.BlockType.Normal_Block;
            this.well_type = Global.WellType.ShutIn;
            this.BHP = new double[steps_memory];
            this.q_oil = new double[steps_memory]; this.q_gas = new double[steps_memory]; this.q_water = new double[steps_memory];
        }

        public void updateProperties(PVT pvt, Kr kr, Porosity porosity, double P, double Sw, double So, double Sg, double p_increment = 0)
        {
            this.dp = Global.epsilon_p;
            this.dso = Global.epsilon_s;
            this.dsw = Global.epsilon_s;
            this.dsg = Global.epsilon_s;

            // pressure
            this.P_previous_step = this.P[0];
            this.P[0] = P; this.P[1] = P + Global.epsilon_p; this.P[2] = this.P[1] + dp;

            // rock
            this.porosity[0] = porosity.getPorosity(P); this.porosity[1] = porosity.getPorosity(this.P[1]); this.porosity[2] = porosity.getPorosity(this.P[2]);

            // fluid
            this.Bo[0] = pvt.getFVF(Global.Phase.Oil, P); this.Bo[1] = pvt.getFVF(Global.Phase.Oil, this.P[1]); this.Bo[2] = pvt.getFVF(Global.Phase.Oil, this.P[2]);
            this.Bw[0] = pvt.getFVF(Global.Phase.Water, P); this.Bw[1] = pvt.getFVF(Global.Phase.Water, this.P[1]); this.Bw[2] = pvt.getFVF(Global.Phase.Water, this.P[2]);
            this.Bg[0] = pvt.getFVF(Global.Phase.Gas, P); this.Bg[1] = pvt.getFVF(Global.Phase.Gas, this.P[1]); this.Bg[2] = pvt.getFVF(Global.Phase.Gas, this.P[2]);

            this.viscosity_oil[0] = pvt.getViscosity(Global.Phase.Oil, P); this.viscosity_oil[1] = pvt.getViscosity(Global.Phase.Oil, this.P[1]); this.viscosity_oil[2] = pvt.getViscosity(Global.Phase.Oil, this.P[2]);
            this.viscosity_water[0] = pvt.getViscosity(Global.Phase.Water, P); this.viscosity_water[1] = pvt.getViscosity(Global.Phase.Water, this.P[1]); this.viscosity_water[2] = pvt.getViscosity(Global.Phase.Water, this.P[2]);
            this.viscosity_gas[0] = pvt.getViscosity(Global.Phase.Gas, P); this.viscosity_gas[1] = pvt.getViscosity(Global.Phase.Gas, this.P[1]); this.viscosity_gas[2] = pvt.getViscosity(Global.Phase.Gas, this.P[2]);

            this.Rso[0] = pvt.getRs(Global.Phase.Oil, P); this.Rso[1] = pvt.getRs(Global.Phase.Oil, this.P[1]); this.Rso[2] = pvt.getRs(Global.Phase.Oil, this.P[2]);

            this.So[0] = So; this.So[1] = So; this.So[2] = this.So[1] + dso;
            this.Sw[0] = Sw; this.Sw[1] = Sw; this.Sw[2] = this.Sw[1] + dsw;
            this.Sg[0] = Sg; this.Sg[1] = Sg; this.Sg[2] = this.Sg[1] + dsg;

            // Kro is only dependent on Sg.
            // if So is at max "Sg = 0" , then So can not increment up.
            this.Kro[0] = kr.getKr(Global.Phase.Oil, Sg); this.Kro[1] = kr.getKr(Global.Phase.Oil, this.Sg[1]); this.Kro[2] = kr.getKr(Global.Phase.Oil, this.Sg[2]);
            // Krw value is always equal to zero in this case.
            this.Krw[0] = 0; this.Krw[1] = 0; this.Krw[2] = 0;
            this.Krg[0] = kr.getKr(Global.Phase.Gas, Sg); this.Krg[1] = kr.getKr(Global.Phase.Gas, this.Sg[1]); this.Krg[2] = kr.getKr(Global.Phase.Gas, this.Sg[2]);

            // volumetric
            this.Vp[0] = this.bulk_volume * this.porosity[0]; this.Vp[1] = this.bulk_volume * this.porosity[1]; this.Vp[2] = this.bulk_volume * this.porosity[2];

            // well blocks
            if (this.type == Global.BlockType.Well_Block)
            {
                this.BHP[0] = Well.WellData.calculatePwf(this, this.P[0], this.Kro[0], this.viscosity_oil[0], this.Bo[0]);
                this.BHP[1] = Well.WellData.calculatePwf(this, this.P[1], this.Kro[1], this.viscosity_oil[1], this.Bo[1]);
                this.BHP[2] = Well.WellData.calculatePwf(this, this.P[2], this.Kro[2], this.viscosity_oil[2], this.Bo[2]);

                this.q_oil[1] = this.q_oil[0];
                this.q_oil[2] = this.q_oil[0];

                //this.q_oil[0] = Well.WellData.calculateFlow_Rate(this.P[0], this.specified_BHP, this.Krg[0], this.viscosity_gas[0], this.WI, this.Bg[0]);
                //this.q_oil[1] = Well.WellData.calculateFlow_Rate(this.P[1], this.specified_BHP, this.Krg[1], this.viscosity_gas[1], this.WI, this.Bg[1]);
                //this.q_oil[2] = Well.WellData.calculateFlow_Rate(this.P[2], this.specified_BHP, this.Krg[2], this.viscosity_gas[2], this.WI, this.Bg[2]);

                this.q_gas[0] = Well.WellData.calculateFlow_Rate(this.P[0], this.BHP[0], this.Krg[0], this.viscosity_gas[0], this.WI, this.Bg[0]);
                this.q_gas[1] = Well.WellData.calculateFlow_Rate(this.P[1], this.BHP[1], this.Krg[1], this.viscosity_gas[1], this.WI, this.Bg[1]);
                this.q_gas[2] = Well.WellData.calculateFlow_Rate(this.P[2], this.specified_BHP, this.Krg[2], this.viscosity_gas[2], this.WI, this.Bg[2]);

                
            }

            total_qg = this.q_gas[0] + this.Rso[0] * this.q_oil[0];
            GOR = total_qg / q_oil[0];

        }

        // this method is not used
        //public void updateProperties_n1(PVT pvt, Kr kr, Porosity porosity, double P, double Sw, double So, double Sg)
        //{
        //    // pressure
        //    this.P[1] = P;

        //    // rock
        //    this.porosity[1] = porosity.getPorosity(P);

        //    // fluid
        //    this.Bo[1] = pvt.getFVF(Global.Phase.Oil, P);
        //    this.Bw[1] = pvt.getFVF(Global.Phase.Water, P);
        //    this.Bg[1] = pvt.getFVF(Global.Phase.Gas, P);

        //    this.viscosity_oil[1] = pvt.getViscosity(Global.Phase.Oil, P);
        //    this.viscosity_water[1] = pvt.getViscosity(Global.Phase.Water, P);
        //    this.viscosity_gas[1] = pvt.getViscosity(Global.Phase.Gas, P);

        //    this.Rso[1] = pvt.getRs(Global.Phase.Oil, P);

        //    this.So[1] = So;
        //    this.Sw[1] = Sw;
        //    this.Sg[1] = Sg;

        //    this.Kro[1] = kr.getKr(Global.Phase.Oil, Sg);
        //    this.Krw[1] = kr.getKr(Global.Phase.Water, Sg);
        //    this.Krg[1] = kr.getKr(Global.Phase.Gas, Sg);

        //    //this.Po = new double[steps_memory]; this.Pg = new double[steps_memory]; this.Pw = new double[steps_memory];

        //    // volumetric
        //    this.Vp[1] = this.bulk_volume * this.porosity[1];

        //    // well
        //    if (this.type == Global.BlockType.Well_Block)
        //    {
        //        this.BHP[1] = Well.WellData.calculatePwf(this, this.P[1], this.Kro[1], this.viscosity_oil[1], this.Bo[1]);

        //        this.q_gas[1] = Well.WellData.calculateFlow_Rate(this.P[1], this.BHP[1], this.Krg[1], this.viscosity_gas[1], this.WI, this.Bg[1]);
        //    }
           
        //}

        public void updateProperties_n1_k1(PVT pvt, Kr kr, Porosity porosity, double P, double Sw, double So, double Sg)
        {
            // time level n+1
            ////////////////////////////////////////////////////////////////

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

            //this.Rso[2] = this.Rso[1];
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
            if (this.type == Global.BlockType.Well_Block)
            {
                this.BHP[1] = Well.WellData.calculatePwf(this, this.P[1], this.Kro[1], this.viscosity_oil[1], this.Bo[1]);
                //Console.WriteLine("BHP : " + BHP[1]);

                this.q_gas[1] = Well.WellData.calculateFlow_Rate(this.P[1], this.BHP[1], this.Krg[1], this.viscosity_gas[1], this.WI, this.Bg[1]);
                //Console.WriteLine("q_gas : " + q_gas[1]);
            }


            // time level k+1
            ////////////////////////////////////////////////////////////////
            //this.dp = Global.epsilon_p * this.P[0];
            //this.dso = Global.epsilon_s * this.So[0];
            //this.dsw = Global.epsilon_s * this.Sw[0];
            ////this.dsg = this.Sg[0] > 0 ? Global.epsilon_s * this.Sg[0] : Global.epsilon_s;
            ////this.dsg = this.dsg < 1E-5 ? 1E-5 : this.dsg;
            //dsg = 1E-4;

            P += dp;

            // pressure
            this.P[2] = P;

            // rock
            this.porosity[2] = porosity.getPorosity(P);

            // fluid
            this.Bo[2] = pvt.getFVF(Global.Phase.Oil, P);
            this.Bw[2] = pvt.getFVF(Global.Phase.Water, P);
            this.Bg[2] = pvt.getFVF(Global.Phase.Gas, P);

            this.viscosity_oil[2] = pvt.getViscosity(Global.Phase.Oil, P);
            this.viscosity_water[2] = pvt.getViscosity(Global.Phase.Water, P);
            this.viscosity_gas[2] = pvt.getViscosity(Global.Phase.Gas, P);

            this.Rso[2] = pvt.getRs(Global.Phase.Oil, P);

            this.Sw[2] = Sw + dsw;
            this.Sg[2] = Sg + dsg;
            this.So[2] = So + dso;

            Sg = this.Sg[2];

            this.Kro[2] = kr.getKr(Global.Phase.Oil, Sg);
            this.Krw[2] = kr.getKr(Global.Phase.Water, Sg);
            this.Krg[2] = kr.getKr(Global.Phase.Gas, Sg);


            // volumetric
            this.Vp[2] = this.bulk_volume * this.porosity[2];

            // well
            if (this.type == Global.BlockType.Well_Block)
            {
                this.BHP[2] = Well.WellData.calculatePwf(this, this.P[2], this.Kro[2], this.viscosity_oil[2], this.Bo[2]);

                this.q_oil[2] = this.q_oil[1];

                this.q_gas[2] = Well.WellData.calculateFlow_Rate(this.P[2], this.BHP[2], this.Krg[2], this.viscosity_gas[2], this.WI, this.Bg[2]);
            }
            
        }

        internal void reset_n1(PVT pvt, Kr kr, Porosity porosity, double p_increment = 0)
        {
            // time level n+1
            ////////////////////////////////////////////////////////////////

            double P = this.P[0] + p_increment;

            // pressure
            this.P[1] = P + Global.epsilon_p;

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

            this.So[1] = this.So[0];
            this.Sw[1] = this.Sw[0];
            this.Sg[1] = this.Sg[0];

            double Sg = this.Sg[1];

            this.Kro[1] = kr.getKr(Global.Phase.Oil, Sg);
            this.Krw[1] = kr.getKr(Global.Phase.Water, Sg);
            this.Krg[1] = kr.getKr(Global.Phase.Gas, Sg);


            // volumetric
            this.Vp[1] = this.bulk_volume * this.porosity[1];

            // well
            if (this.type == Global.BlockType.Well_Block)
            {
                this.BHP[1] = Well.WellData.calculatePwf(this, this.P[1], this.Kro[1], this.viscosity_oil[1], this.Bo[1]);

                this.q_gas[1] = Well.WellData.calculateFlow_Rate(this.P[1], this.BHP[1], this.Krg[1], this.viscosity_gas[1], this.WI, this.Bg[1]);
            }


            // time level k+1
            ////////////////////////////////////////////////////////////////
            //this.dp = Global.epsilon_p * this.P[0];
            //this.dso = Global.epsilon_s * this.So[0];
            //this.dsw = Global.epsilon_s * this.Sw[0];
            ////this.dsg = this.Sg[0] > 0 ? Global.epsilon_s * this.Sg[0] : Global.epsilon_s;
            ////this.dsg = this.dsg < 1E-5 ? 1E-5 : this.dsg;
            //dsg = 1E-4;

            P = this.P[1] + dp;

            // pressure
            this.P[2] = P;

            // rock
            this.porosity[2] = porosity.getPorosity(P);

            // fluid
            this.Bo[2] = pvt.getFVF(Global.Phase.Oil, P);
            this.Bw[2] = pvt.getFVF(Global.Phase.Water, P);
            this.Bg[2] = pvt.getFVF(Global.Phase.Gas, P);

            this.viscosity_oil[2] = pvt.getViscosity(Global.Phase.Oil, P);
            this.viscosity_water[2] = pvt.getViscosity(Global.Phase.Water, P);
            this.viscosity_gas[2] = pvt.getViscosity(Global.Phase.Gas, P);

            this.Rso[2] = pvt.getRs(Global.Phase.Oil, P);

            this.So[2] = this.So[1] + dso;
            this.Sw[2] = this.Sw[1] + dsw;
            this.Sg[2] = this.Sg[1] + dsg;

            Sg = this.Sg[2];

            this.Kro[2] = kr.getKr(Global.Phase.Oil, Sg);
            this.Krw[2] = kr.getKr(Global.Phase.Water, Sg);
            this.Krg[2] = kr.getKr(Global.Phase.Gas, Sg);


            // volumetric
            this.Vp[2] = this.bulk_volume * this.porosity[2];

            // well
            if (this.type == Global.BlockType.Well_Block)
            {
                this.BHP[2] = Well.WellData.calculatePwf(this, this.P[2], this.Kro[2], this.viscosity_oil[2], this.Bo[2]);

                this.q_gas[2] = Well.WellData.calculateFlow_Rate(this.P[2], this.BHP[2], this.Krg[2], this.viscosity_gas[2], this.WI, this.Bg[2]);
            }
            
        }

        public double calculateR(BaseBlock[] grid, Global.Phase phase, double time_step, bool solubleGasPresent = true)
        {
            // Note : Kro, Krw, Krg are all calculated based on Sg.

            BaseBlock block = this;

            BaseBlock upstream_block, downstream_block, neighbour_block;
            double transmissibility;

            double Kr, B, viscosity, Rso;
            double R = 0;

            double temp, accumulation_term, production_term = 0;


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

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                production_term = block.q_oil[1];
                block.production_term_oil = production_term;

                R -= (accumulation_term + production_term);
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

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                production_term = block.q_water[1];
                block.production_term_water = production_term;


                R -= (accumulation_term + production_term);
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

                    if (solubleGasPresent)
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

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    accumulation_term += 1 / (Global.a * time_step) * ((block.Rso[1] * block.Vp[1] * (1 - block.Sw[1] - block.Sg[1]) / block.Bo[1]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                block.accumulation_term_gas = accumulation_term;


                if (block.type == Global.BlockType.Well_Block)
                {
                    if (block.well_type == Global.WellType.Production)
                    {
                        production_term = block.q_gas[1];
                        // check for presence of soluble_gas in simulation_data
                        if (solubleGasPresent)
                        {
                            production_term += block.Rso[1] * block.q_oil[1] ;
                        }
                    }
                    else if (block.well_type == Global.WellType.Injection)
                    {
                        production_term = block.specified_flow_rate;
                    }
                }
                block.production_term_gas = production_term;


                R -= (accumulation_term + production_term);
            }

            return R;

        }

        public double calculateDerivative(BaseBlock[] grid, Global.Phase phase, double time_step, int neighbour_block_index, Global.Variable variable, bool solubleGasPresent = true)
        {
            BaseBlock block = this;

            BaseBlock upstream_block, downstream_block, neighbour_block;
            double transmissibility;

            double Kr, B, viscosity, Rso;
            double R = 0;

            double temp, accumulation_term, production_term = 0;
            int s = 1, ns = 1;


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
                    ns = 2;
                }
                else
                {
                    sgd = 2;
                    s = 2;
                }
            }
            else
            {
                if (!this_block)
                {
                    nswd = 2;
                    ns = 2;
                }
                else
                {
                    swd = 2;
                    s = 2;
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

                        if (block.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;

                            Kr = upstream_block.Kro[s];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;

                            Kr = upstream_block.Kro[ns];
                        }

                        B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                        R += temp;

                    }

                    accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * (1 - block.Sw[swd] - block.Sg[sgd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                    production_term = block.production_term_oil;

                    R -= (accumulation_term + production_term);

                }
                #endregion
                #region NextBlock
                else
                {
                    neighbour_block = grid[block.neighbour_blocks_indices[neighbour_block_index]];
                    transmissibility = block.transmissibility_list[neighbour_block_index];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;

                        Kr = upstream_block.Kro[s];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;

                        Kr = upstream_block.Kro[ns];
                    }

                    B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                    viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                    R = temp;

                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i != neighbour_block_index)
                        {
                            R += block.transmissibility_terms_oil[i];
                        }
                    }

                    R -= (block.accumulation_term_oil + block.production_term_oil);
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

                        if (block.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;

                            Kr = upstream_block.Krw[sgd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;

                            Kr = upstream_block.Krw[nsgd];
                        }

                        B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                        viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                        R += temp;
                    }

                    accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * (block.Sw[swd]) / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0])) + block.q_water[1];

                    R -= accumulation_term;
                }
                #endregion
                #region NextBlock
                else
                {
                    neighbour_block = grid[block.neighbour_blocks_indices[neighbour_block_index]];
                    transmissibility = block.transmissibility_list[neighbour_block_index];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;

                        Kr = upstream_block.Krw[sgd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;

                        Kr = upstream_block.Krw[nsgd];
                    }

                    B = 0.5 * (block.Bw[pd] + neighbour_block.Bw[npd]);
                    viscosity = 0.5 * (block.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                    R += temp;

                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i != neighbour_block_index)
                        {
                            R += block.transmissibility_terms_water[i];
                        }
                    }

                    R -= block.accumulation_term_water;
                    //R -= block.production_term_water;
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

                        if (block.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = block;
                            downstream_block = neighbour_block;

                            Kr = upstream_block.Krg[sgd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = block;

                            Kr = upstream_block.Krg[nsgd];
                        }

                        B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                        viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);

                        temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                        if (solubleGasPresent)
                        {
                            if (upstream_block.index == block.index)
                            {
                                Kr = upstream_block.Kro[sgd];
                            }
                            else
                            {
                                Kr = upstream_block.Kro[nsgd];
                            }

                            B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                            viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);
                            Rso = 0.5 * (block.Rso[pd] + neighbour_block.Rso[npd]);

                            temp += Rso * transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);
                        }

                        R += temp;
                    }

                    accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * block.Sg[sgd] / block.Bg[pd]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                    // check for presence of soluble_gas in simulation_data
                    if (solubleGasPresent)
                    {
                        accumulation_term += 1 / (Global.a * time_step) * ((block.Rso[pd] * block.Vp[pd] * (1 - block.Sw[swd] - block.Sg[sgd]) / block.Bo[pd]) - (block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0]));
                    }

                    if (block.type == Global.BlockType.Well_Block)
                    {
                        if (block.well_type == Global.WellType.Production)
                        {
                            double bhp = Well.WellData.calculatePwf(block, block.P[pd], block.Kro[sgd], block.viscosity_oil[pd], block.Bo[pd]);
                            production_term += Well.WellData.calculateFlow_Rate(block.P[pd], bhp, block.Krg[sgd], block.viscosity_gas[pd], block.WI, block.Bg[pd]);
                            // check for presence of soluble_gas in simulation_data
                            if (solubleGasPresent)
                            {
                                production_term += block.Rso[pd] * block.q_oil[1];
                            }
                        }
                        else if (block.well_type == Global.WellType.Injection)
                        {
                            //To-Do: implement gas injection treatment.
                            //block.q_gas[2] = block.q_gas[1];
                            production_term = block.specified_flow_rate;
                        }
                    }


                    R -= (accumulation_term + production_term);
                }
                #endregion
                #region NextBlock
                else
                {
                    neighbour_block = grid[block.neighbour_blocks_indices[neighbour_block_index]];
                    transmissibility = block.transmissibility_list[neighbour_block_index];

                    if (block.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = block;
                        downstream_block = neighbour_block;

                        Kr = upstream_block.Krg[sgd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = block;

                        Kr = upstream_block.Krg[nsgd];
                    }

                    B = 0.5 * (block.Bg[pd] + neighbour_block.Bg[npd]);
                    viscosity = 0.5 * (block.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);

                    temp = transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);

                    if (solubleGasPresent)
                    {
                        if (upstream_block.index == block.index)
                        {
                            Kr = upstream_block.Kro[sgd];
                        }
                        else
                        {
                            Kr = upstream_block.Kro[nsgd];
                        }

                        B = 0.5 * (block.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (block.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);
                        Rso = 0.5 * (block.Rso[pd] + neighbour_block.Rso[npd]);

                        temp += Rso * transmissibility * Kr / (viscosity * B) * (neighbour_block.P[npd] - block.P[pd]);
                    }

                    R += temp;

                    for (int i = 0; i < block.neighbour_blocks_indices.Length; i++)
                    {
                        if (i != neighbour_block_index)
                        {
                            R += block.transmissibility_terms_gas[i];
                        }
                    }

                    R -= (block.accumulation_term_gas + block.production_term_gas);
                }
                #endregion
            }
            #endregion

            return R;

        }

        public double calculateJacobian_New(SimulationData data, Global.Phase equation_phase, int neighbour_block_index, Global.Variable variable)
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


            if (equation_phase == Global.Phase.Oil)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < this.neighbour_blocks_indices.Length; i++)
                    {
                        if (this.neighbour_blocks_indices[i] == -1)
                        {
                            continue;
                        }

                        neighbour_block = data.grid[this.neighbour_blocks_indices[i]];

                        if (this.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = this;
                            downstream_block = neighbour_block;

                            kr = upstream_block.Kro[sgd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = this;

                            kr = upstream_block.Kro[nsgd];
                        }

                        B = 0.5 * (this.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (this.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);
                        transmissiblity = this.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);

                        //this.transmissibility_terms_oil[i] = transmissibility_temp;

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.time_step * Global.a) * ((this.Vp[pd] * (1 - this.Sg[sgd] - this.Sw[swd]) / this.Bo[pd]) - (this.Vp[0] * this.So[0] / this.Bo[0]));

                    //this.accumulation_term_oil = accumulation_term;
                    // production term
                    production_term = this.q_oil[1];

                    //this.production_term_oil = production_term;
                }
                else
                {
                    // transmissibility term
                    neighbour_block = data.grid[this.neighbour_blocks_indices[neighbour_block_index]];

                    if (this.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = this;
                        downstream_block = neighbour_block;

                        kr = upstream_block.Kro[sgd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = this;

                        kr = upstream_block.Kro[nsgd];
                    }

                    B = 0.5 * (this.Bo[pd] + neighbour_block.Bo[npd]);
                    viscosity = 0.5 * (this.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);
                    transmissiblity = this.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < this.neighbour_blocks_indices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += this.transmissibility_terms_oil[i];
                    }
                    // accumulation term
                    accumulation_term = this.accumulation_term_oil;
                    // production term
                    production_term = this.production_term_oil;

                }
            }
            else if (equation_phase == Global.Phase.Gas)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < this.neighbour_blocks_indices.Length; i++)
                    {
                        if (this.neighbour_blocks_indices[i] == -1)
                        {
                            continue;
                        }

                        neighbour_block = data.grid[this.neighbour_blocks_indices[i]];

                        if (this.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = this;
                            downstream_block = neighbour_block;

                            kr = upstream_block.Krg[sgd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = this;

                            kr = upstream_block.Krg[nsgd];
                        }

                        B = 0.5 * (this.Bg[pd] + neighbour_block.Bg[npd]);
                        viscosity = 0.5 * (this.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);
                        transmissiblity = this.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);

                        if (data.solubleGasPresent)
                        {
                            if (upstream_block.index == this.index)
                            {
                                kr = upstream_block.Kro[sgd];
                            }
                            else
                            {
                                kr = upstream_block.Kro[nsgd];
                            }

                            Rso = 0.5 * (this.Rso[pd] + neighbour_block.Rso[npd]);
                            B = 0.5 * (this.Bo[pd] + neighbour_block.Bo[npd]);
                            viscosity = 0.5 * (this.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                            transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);
                        }

                        //this.transmissibility_terms_gas[i] = transmissibility_temp;

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.time_step * Global.a) * ((this.Vp[pd] * this.Sg[sgd] / this.Bg[pd]) - (this.Vp[0] * this.Sg[0] / this.Bg[0]));

                    if (data.solubleGasPresent)
                    {
                        accumulation_term += 1 / (data.time_step * Global.a) * ( this.Rso[pd] * (this.Vp[pd] * (1 - this.Sg[sgd] - this.Sw[swd]) / this.Bo[pd]) - this.Rso[0] * (this.Vp[0] * this.So[0] / this.Bo[0]));
                    }

                    //this.accumulation_term_gas = accumulation_term;
                    // production term
                    if (this.type == Global.BlockType.Well_Block)
                    {
                        if (this.well_type == Global.WellType.Production)
                        {
                            production_term = this.q_gas[1];

                            if (data.solubleGasPresent)
                            {
                                production_term += this.Rso[pd] * this.q_oil[1];
                            }
                        }
                        else
                        {
                            production_term = this.specified_flow_rate;
                        }
                    }
                    

                    //this.production_term_gas = production_term;
                }
                else
                {
                    // transmissibility term
                    neighbour_block = data.grid[this.neighbour_blocks_indices[neighbour_block_index]];

                    if (this.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = this;
                        downstream_block = neighbour_block;

                        kr = upstream_block.Krg[sgd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = this;

                        kr = upstream_block.Krg[nsgd];
                    }

                    B = 0.5 * (this.Bg[pd] + neighbour_block.Bg[npd]);
                    viscosity = 0.5 * (this.viscosity_gas[pd] + neighbour_block.viscosity_gas[npd]);
                    transmissiblity = this.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);

                    if (data.solubleGasPresent)
                    {
                        if (upstream_block.index == this.index)
                        {
                            kr = upstream_block.Kro[sgd];
                        }
                        else
                        {
                            kr = upstream_block.Kro[nsgd];
                        }

                        Rso = 0.5 * (this.Rso[pd] + neighbour_block.Rso[npd]);
                        B = 0.5 * (this.Bo[pd] + neighbour_block.Bo[npd]);
                        viscosity = 0.5 * (this.viscosity_oil[pd] + neighbour_block.viscosity_oil[npd]);

                        transmissibility_temp += Rso * transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);
                    }

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < this.neighbour_blocks_indices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += this.transmissibility_terms_gas[i];
                    }
                    // accumulation term
                    accumulation_term = this.accumulation_term_gas;
                    // production term
                    production_term = this.production_term_gas;
                }
            }
            else if (equation_phase == Global.Phase.Water)
            {
                if (this_block)
                {
                    // transmissibility term
                    for (int i = 0; i < this.neighbour_blocks_indices.Length; i++)
                    {
                        if (this.neighbour_blocks_indices[i] == -1)
                        {
                            continue;
                        }

                        neighbour_block = data.grid[this.neighbour_blocks_indices[i]];

                        if (this.P[1] >= neighbour_block.P[1])
                        {
                            upstream_block = this;
                            downstream_block = neighbour_block;

                            kr = upstream_block.Krw[swd];
                        }
                        else
                        {
                            upstream_block = neighbour_block;
                            downstream_block = this;

                            kr = upstream_block.Krw[nswd];
                        }

                        B = 0.5 * (this.Bw[pd] + neighbour_block.Bw[npd]);
                        viscosity = 0.5 * (this.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);
                        transmissiblity = this.transmissibility_list[i];

                        transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);

                        //this.transmissibility_terms_water[i] = transmissibility_temp;

                        transmissibility_term += transmissibility_temp;
                    }
                    // accumulation term
                    accumulation_term = 1 / (data.time_step * Global.a) * ((this.Vp[pd] * this.Sw[swd] / this.Bw[pd]) - (this.Vp[0] * this.Sw[0] / this.Bw[0]));

                    //this.accumulation_term_water = accumulation_term;
                    // production term
                    production_term = this.q_water[1];

                    //this.production_term_water = production_term;
                }
                else
                {
                    // transmissibility term
                    neighbour_block = data.grid[this.neighbour_blocks_indices[neighbour_block_index]];

                    if (this.P[1] >= neighbour_block.P[1])
                    {
                        upstream_block = this;
                        downstream_block = neighbour_block;

                        kr = upstream_block.Krw[swd];
                    }
                    else
                    {
                        upstream_block = neighbour_block;
                        downstream_block = this;

                        kr = upstream_block.Krw[nswd];
                    }

                    B = 0.5 * (this.Bw[pd] + neighbour_block.Bw[npd]);
                    viscosity = 0.5 * (this.viscosity_water[pd] + neighbour_block.viscosity_water[npd]);
                    transmissiblity = this.transmissibility_list[neighbour_block_index];

                    transmissibility_temp = transmissiblity * kr / (viscosity * B) * (neighbour_block.P[npd] - this.P[pd]);

                    transmissibility_term = transmissibility_temp;
                    for (int i = 0; i < this.neighbour_blocks_indices.Length; i++)
                    {
                        if (i == neighbour_block_index)
                        {
                            continue;
                        }

                        transmissibility_term += this.transmissibility_terms_water[i];
                    }
                    // accumulation term
                    accumulation_term = this.accumulation_term_water;
                    // production term
                    production_term = this.production_term_water;
                }
            }


            R_plus = transmissibility_term - accumulation_term - production_term;
            return R_plus;
        }

        public double single_block_R(SimulationData data, Global.Phase phase)
        {
            // Note : Kro, Krw, Krg are all calculated based on Sg.

            BaseBlock block = this;

            BaseBlock upstream_block, downstream_block, neighbour_block;
            double transmissibility;

            double Kr, B, viscosity, Rso;
            double R = 0;

            double temp, accumulation_term, production_term = 0;

            bool solubleGasPresent = data.solubleGasPresent;
            double time_step = data.time_step;


            if (phase == Global.Phase.Oil)
            {
                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * (1 - block.Sg[1]) / block.Bo[1]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                production_term = block.q_oil[1];
                //production_term = Well.WellData.calculateFlow_Rate(this.P[1], this.specified_BHP, this.Kro[1], this.viscosity_oil[1], this.WI, this.Bo[1]);
                block.production_term_oil = production_term;

                R -= (accumulation_term + production_term);
            }
            else if (phase == Global.Phase.Water)
            {
                

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.Sw[1] / block.Bw[1]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                R -= (accumulation_term + production_term);
            }
            else
            {

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[1] * block.Sg[1] / block.Bg[1]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    accumulation_term += 1 / (Global.a * time_step) * ((block.Rso[1] * data.pvt.Rs_SG(this.Sg[1]) * block.Vp[1] * (1 - block.Sg[1]) / block.Bo[1]) - (block.Rso[0] * data.pvt.Rs_SG(this.Sg[0]) * block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                block.accumulation_term_gas = accumulation_term;


                production_term = block.q_gas[1];
                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    production_term += block.Rso[1] * data.pvt.Rs_SG(this.Sg[1]) * block.production_term_oil;
                }
                block.production_term_gas = production_term;


                R -= (accumulation_term + production_term);
            }

            return R;
        }

        public double single_block_perturb(SimulationData data, Global.Phase phase, int neighbour_block_index, Global.Variable variable)
        {
            BaseBlock block = this;
            double time_step = data.time_step;
            bool solubleGasPresent = data.solubleGasPresent;
            double R = 0;

            double R_plus = 0;

            double transmissibility_temp = 0;

            double transmissibility_term = 0;
            double accumulation_term = 0;
            double production_term = 0;

            double kr, viscosity, B, Rso;
            BaseBlock upstream_block, downstream_block, neighbour_block;

            double transmissiblity;

            double COP = 0; double COG = 0; double CGP = 0; double CGG = 0;

            #region Variable Dependencies
            int pd = 1, swd = 1, sgd = 1;


            if (variable == Global.Variable.Pressure)
            {
                pd = 2;
            }
            else if (variable == Global.Variable.Saturation_Gas)
            {
                sgd = 2;
            }
            else if (variable == Global.Variable.Saturation_Water)
            {
                swd = 2;
            }
            #endregion

            if (phase == Global.Phase.Oil)
            {
                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * (1 - block.Sg[sgd]) / block.Bo[pd]) - (block.Vp[0] * block.So[0] / block.Bo[0]));
                block.accumulation_term_oil = accumulation_term;

                COP = 1 / (Global.a * time_step) * (this.phi_dash() / this.Bo[0] + this.Vp[1] * this.FVF_dash(Global.Phase.Oil)) * (1 - this.Sg[0]);
                COG = 1 / (Global.a * time_step)  * - 1 * this.Vp[1] / this.Bo[1];

                //accumulation_term = 1 / (Global.a * time_step) * (COP * (this.P[pd] - this.P[0]) + COG * (this.Sg[sgd] - this.Sg[0]));

                production_term = block.q_oil[1];
                block.production_term_oil = production_term;

                R -= (accumulation_term + production_term);
            }
            else if (phase == Global.Phase.Water)
            {


                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * block.Sw[swd] / block.Bw[pd]) - (block.Vp[0] * block.Sw[0] / block.Bw[0]));
                block.accumulation_term_water = accumulation_term;

                R -= (accumulation_term + production_term);
            }
            else
            {

                accumulation_term = 1 / (Global.a * time_step) * ((block.Vp[pd] * block.Sg[sgd] / block.Bg[pd]) - (block.Vp[0] * block.Sg[0] / block.Bg[0]));
                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    accumulation_term += 1 / (Global.a * time_step) * ((block.Rso[pd] * data.pvt.Rs_SG(this.Sg[sgd]) * block.Vp[pd] * (1 - block.Sg[sgd]) / block.Bo[pd]) - (block.Rso[0] * data.pvt.Rs_SG(this.Sg[0]) * block.Vp[0] * block.So[0] / block.Bo[0]));
                }
                block.accumulation_term_gas = accumulation_term;

                CGP = 1 / (Global.a * time_step) * (((this.phi_dash() / this.Bo[0] + this.Vp[1] * this.FVF_dash(Global.Phase.Oil)) * this.Rso[0] + (this.Vp[1] / this.Bo[1]) * this.Rso_dash()) * (1 - this.Sg[0]) + (this.phi_dash() / this.Bg[0] + this.Vp[1] * this.FVF_dash(Global.Phase.Gas)) * this.Sg[0]);
                CGG = 1 / (Global.a * time_step) * (this.Vp[1] / this.Bg[1] - this.Vp[1] / this.Bo[1] * this.Rso[1]);

                //accumulation_term = 1 / (Global.a * time_step) * (CGP * (this.P[pd] - this.P[0]) + CGG * (this.Sg[sgd] - this.Sg[0]));

                production_term = block.q_gas[1];
                // check for presence of soluble_gas in simulation_data
                if (solubleGasPresent)
                {
                    production_term += block.Rso[pd] * data.pvt.Rs_SG(this.Sg[sgd]) * block.q_oil[1];
                }
                block.production_term_gas = production_term;


                R -= (accumulation_term + production_term);
            }

            return R;
        }

        public double qasim_jacobian(SimulationData data, Global.Phase phase, int neighbour_block_index, Global.Variable variable)
        {
            double temp = 0;

            double COP = 0, COG = 0, CGP = 0, CGG = 0, q = 0;

            if (phase == Global.Phase.Oil)
            {
                if (variable == Global.Variable.Pressure)
                {
                    double phi_dash = this.phi_dash();
                    double Bo_dash = this.FVF_dash(Global.Phase.Oil);

                    COP = 1 / (Global.a * data.time_step) * (this.phi_dash() / this.Bo[0] + this.Vp[1] * this.FVF_dash(Global.Phase.Oil)) * (1 - this.Sg[0]);
                    return -COP;
                }
                else
                {
                    COG = -1 / (Global.a * data.time_step) * (this.Vp[1] / this.Bo[1]);
                    return -COG;
                }
            }
            if (phase == Global.Phase.Gas)
            {
                if (variable == Global.Variable.Pressure)
                {
                    CGP = 1 / (Global.a * data.time_step) * (((this.phi_dash() / this.Bo[0] + this.Vp[1] * this.FVF_dash(Global.Phase.Oil)) * this.Rso[0] + (this.Vp[1] / this.Bo[1]) * this.Rso_dash()) * (1 - this.Sg[0]) + (this.phi_dash() / this.Bg[0] + this.Vp[1] * this.FVF_dash(Global.Phase.Gas)) * this.Sg[0]);
                    q = (this.Rso[1] - this.Rso[2]) / (this.P[1] - this.P[2]) * this.q_oil[1];

                    return -CGP - q;
                }
                else
                {
                    CGG = 1 / (Global.a * data.time_step) * (this.Vp[1] / this.Bg[1] - (this.Vp[1] / this.Bo[1]) * this.Rso[1]);
                    return -CGG;
                }
            }
            return temp;
        }

        public double phi_dash()
        {
            double P_difference = this.P[1] - this.P[0];
            return (this.Vp[1] - this.Vp[0]) / (P_difference);
        }

        public double Rso_dash()
        {
            double P_difference = this.P[1] - this.P[0];
            return (this.Rso[1] - this.Rso[0]) / (P_difference);
        }

        public double FVF_dash(Global.Phase phase)
        {
            double P_difference = this.P[1] - this.P[0];

            double FVF_1 = 1, FVF_0 = 1;
            if (phase == Global.Phase.Oil)
            {
                FVF_1 = this.Bo[1];
                FVF_0 = this.Bo[0];
            }
            else if (phase == Global.Phase.Water)
            {
                FVF_1 = this.Bw[1];
                FVF_0 = this.Bw[0];
            }
            else
            {
                FVF_1 = this.Bg[1];
                FVF_0 = this.Bg[0];
            }
            return (1 / FVF_1 - 1 / FVF_0) / (P_difference);
        }

    }
}