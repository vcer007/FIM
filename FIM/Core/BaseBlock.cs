using System;

/// <summary>
/// The base namespace for the simulator.
/// </summary>
namespace FIM.Core
{
    /// <summary>
    /// The base class that contains the minimum set of common properties that any grid block should have.
    /// </summary>
    /// <remarks>
    ///  <para>This class is the common foundation for storing data of the petro-physical properties of each block 
    ///  regardless of specific block type.</para>
    ///  <para>In order to extend the functionality of the simulator, extension methods should be used 
    ///  to keep the clean structure of the <c>BaseBlock</c> intact.</para>
    ///  <see cref="FIM.Extensions"/>
    ///  <para>To add more properties, consider creating a new class inheriting from <c>BaseBlock</c> 
    ///  to keep the original code at its minimal state 
    ///  and to prevent unnecessary memory storage increase in cases where using the additional properties is not needed.</para>
    /// </remarks>
    public class BaseBlock
    {
        #region Indexing

        /// <summary>
        /// The index of the block in the <see cref="SimulationData.grid"/> array.
        /// </summary>
        public int index;

        /// <summary>
        /// For a structured grid, this property is used to store the layer of the block.
        /// </summary>
        public int layer;

        /// <summary>
        /// A list of the indices of neighbouring blocks.
        /// </summary>
        /// <remarks>
        /// <para>This <see cref="Array"/> contains the indices of the blocks neighboring this block.</para>
        /// <para>The list begins with the index of the top block, then the bottom block then all the other blocks.</para>
        /// <para>A value of -1 is used in cases where the top or bottom blocks are non-present.</para>
        /// </remarks>
        public int[] neighbour_blocks_indices;

        #endregion

        #region data.porosity properties

        /// <summary>
        /// An <see cref="Array"/> of the grid block porosities at different time levels.
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] porosity;

        /// <summary>
        /// An <see cref="Array"/> of the grid block permeabilities in the different directions.
        /// </summary>
        public double[] permeability_list;

        #endregion

        #region Fluid properties

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid saturations at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] Sw, So, Sg;

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid relative permeabilities at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] Krw, Kro, Krg;

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid FVF "Formation Volume Factor" at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] Bw, Bo, Bg;

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid viscosities at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] viscosity_water, viscosity_oil, viscosity_gas;

        /// <summary>
        /// <para>An <see cref="Array"/> of the solution Gas/Oil ratio at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] Rso;

        /// <summary>
        /// <para>An <see cref="Array"/> of the pressure ratio at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.steps_memory"/>
        public double[] P;

        // phase specific properties. used when implementing capillary pressure.
        //public double[] Pw, Po, Pg;

        /// <summary>
        /// The pressure of the previous time step "n-1 time level".
        /// </summary>
        public double P_previous_step;

        #endregion

        #region Volumetric data

        public double[] area_list, delta_x_list, boundary_length_list;

        public double h, bulk_volume;

        public double[] Vp;

        #endregion

        #region Well data
        public Global.BlockType type;
        //public Global.WellType well_type;

        //public double r_equivalent, WI;

        //public double specified_flow_rate, specified_BHP;

        //public double[] BHP;

        //public double well_radius, skin;

        //public double[] q_water, q_oil, q_gas;

        public double total_qg, GOR;

        #endregion

        #region Miscellaneous

        // list of the blocks transmissibilities, brginning with the top then bottom ones.

        public double[] transmissibility_list;


        // these variables are used to save repeated calculations unnecessarily.
        public double[] transmissibility_terms_oil, transmissibility_terms_water, transmissibility_terms_gas;
        public double accumulation_term_oil, accumulation_term_water, accumulation_term_gas;
        public double production_term_oil, production_term_water, production_term_gas;

        public double dp, dso, dsw, dsg;

        #endregion


        /// <summary>
        /// Initializes a new instance of the <see cref="BaseBlock"/> class.
        /// </summary>
        public BaseBlock()
        {
            int steps_memory = Global.steps_memory;

            // data.porosity_calculator
            this.porosity = new double[steps_memory];

            // fluid
            this.Bo = new double[steps_memory]; this.Bg = new double[steps_memory]; this.Bw = new double[steps_memory];

            this.So = new double[steps_memory]; this.Sg = new double[steps_memory]; this.Sw = new double[steps_memory];

            this.Kro = new double[steps_memory]; this.Krg = new double[steps_memory]; this.Krw = new double[steps_memory];

            this.viscosity_oil = new double[steps_memory]; this.viscosity_gas = new double[steps_memory]; this.viscosity_water = new double[steps_memory];

            this.Rso = new double[steps_memory + 1];

            this.P = new double[steps_memory + 1];

            // volumetric
            this.Vp = new double[steps_memory];

            // well
            this.type = Global.BlockType.Normal_Block;

            //this.well_type = Global.WellType.ShutIn;

            //this.BHP = new double[steps_memory];

            //this.q_oil = new double[steps_memory]; this.q_gas = new double[steps_memory]; this.q_water = new double[steps_memory];
        }

        /// <summary>
        /// Updates the properties of the <see cref="BaseBlock"/>.
        /// </summary>
        /// <remarks>
        /// <para>Different time levels properties are updated recursively down-up.</para>
        /// <para>If time level is equal to zero, then n0, n1 and n2, if available, time levels will be updated.</para>
        /// <para>If time level is equal to 1, then only n1 and n2, if available, time levels will be updated.</para>
        /// </remarks>
        /// <param name="data">The <seealso cref="SimulationData"/> object thaat contains all the input data.</param>
        /// <param name="P">The new pressure value.</param>
        /// <param name="Sw">The new water saturation value.</param>
        /// <param name="Sg">The new gas saturation value.</param>
        /// <param name="time_level"><para>Default value is 0.</para>
        /// <para>The time_level at which the properties should be updated.</para>
        /// </param>
        /// <seealso cref="Global.steps_memory"/>
        public void updateProperties(SimulationData data, double P, double Sw, double Sg, int time_level = 0)
        {
            // pressure
            this.P_previous_step = this.P[0];

            this.P[time_level] = P;

            // data.porosity_calculator
            this.porosity[time_level] = data.porosity_calculator.getPorosity(P);

            // fluid
            this.Bo[time_level] = data.pvt.getFVF(Global.Phase.Oil, P);
            this.Bw[time_level] = data.pvt.getFVF(Global.Phase.Water, P);
            this.Bg[time_level] = data.pvt.getFVF(Global.Phase.Gas, P);

            this.viscosity_oil[time_level] = data.pvt.getViscosity(Global.Phase.Oil, P);
            this.viscosity_water[time_level] = data.pvt.getViscosity(Global.Phase.Water, P);
            this.viscosity_gas[time_level] = data.pvt.getViscosity(Global.Phase.Gas, P);

            this.Rso[time_level] = data.pvt.getRs(Global.Phase.Oil, P);

            this.So[time_level] = 1 - Sw - Sg;
            this.Sw[time_level] = Sw;
            this.Sg[time_level] = Sg;

            // Kro is only dependent on Sg.
            this.Kro[time_level] = data.scal.getKr(Global.Phase.Oil, Sg);
            this.Krw[time_level] = data.scal.getKr(Global.Phase.Water, Sg);
            this.Krg[time_level] = data.scal.getKr(Global.Phase.Gas, Sg);

            // volumetric
            this.Vp[time_level] = this.bulk_volume * this.porosity[time_level];

            // well blocks
            //if (this.type == Global.BlockType.Well_Block)
            //{
            //    this.BHP[time_level] = Well.WellData.calculatePwf(this, this.P[time_level], this.Kro[time_level], this.viscosity_oil[time_level], this.Bo[time_level]);

            //    this.q_oil[1] = this.q_oil[0];

            //    this.q_gas[time_level] = Well.WellData.calculateFlow_Rate(this.P[time_level], this.BHP[time_level], this.Krg[time_level], this.viscosity_gas[time_level], this.WI, this.Bg[time_level]);
            //}

            if (type == Global.BlockType.Well_Block)
            {
                for (int i = 0; i < data.wells.Length; i++)
                {
                    if (data.wells[i].index == index)
                    {
                        data.wells[i].update(time_level, this);
                    }
                }
            }

            if (time_level == 0)
            {
                //total_qg = this.q_gas[time_level] + this.Rso[time_level] * this.q_oil[time_level];
                //GOR = total_qg / q_oil[time_level];

                // this way we update n1 time level properties also automatically
                // after updating n0 time level properties.
                // note that the initial guess here for n1 is assumed to be the same values as the ones at n0.
                updateProperties(data, P, Sw, Sg, 1);
            }
            else if (time_level == 1)
            {
                // check if numerical perturbation method is used.
                // this is indicated by three time levels storage "where the third one, n2, is used solely for perturbation".
                if (Global.steps_memory == 3)
                {
                    updateProperties(data, P + Global.epsilon, Sw + Global.epsilon, Sg + Global.epsilon, 2);
                }
            }
        }

        /// <summary>
        /// Resets n1, and n2 time level if available, properties to n0.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <seealso cref="SimulationData"/>
        /// <seealso cref="SimulationData.time_step_slashing_factor"/>
        public void reset_n1(SimulationData data)
        {
            // to reset the properties of the block, simply call updateProperties_n0_n1 with input parameters from
            // the n0 time level.
            updateProperties(data, this.P[0], this.Sw[0], this.Sg[0], 0);
        }

        
    }
}