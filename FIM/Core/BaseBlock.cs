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
        public int[] neighborBlocksIndices;

        #endregion

        #region Rock properties

        /// <summary>
        /// An <see cref="Array"/> of the grid block porosities at different time levels.
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] porosity;

        /// <summary>
        /// An <see cref="Array"/> of the grid block permeabilities in the different directions.
        /// </summary>
        public double[] permeability;

        #endregion

        #region Fluid properties

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid saturations at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] Sw, So, Sg;

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid relative permeabilities at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] Krw, Kro, Krg;

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid FVF "Formation Volume Factor" at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] Bw, Bo, Bg;

        /// <summary>
        /// <para>An <see cref="Array"/> of the fluid viscosities at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] viscosityWater, viscosityOil, viscosityGas;

        /// <summary>
        /// <para>An <see cref="Array"/> of the solution Gas/Oil ratio at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] Rso;

        /// <summary>
        /// <para>An <see cref="Array"/> of the pressure ratio at different time levels.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] P;

        /// <summary>
        /// <para>An <see cref="Array"/> of the capillary pressures for oil/water and gas/oil.</para>
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] Pcow, Pcgo;

        /// <summary>
        /// The pressure of the previous time step "n-1 time level".
        /// </summary>
        public double P_previousStep;

        #endregion

        #region Volumetric data
        /// <summary>
        /// An <see cref="Array"/> of the top, bottom and side faces areas of the block "in this order".
        /// </summary>
        /// <remarks>
        /// The size the array is equal to 2 + number of sideways neighboring blocks.
        /// This is due to the fact that the PeBi grid system implemented depends on flat blocks.
        /// So, each block will have a single top and a single bottom block except if it is a top/bottom boundary block.
        /// </remarks>
        public double[] areaList;

        /// <summary>
        /// An <see cref="Array"/> of the distances from the center of the block to the top, bottom and side faces of the block "in this order".
        /// </summary>
        /// <remarks>
        /// The size the array is equal to 2 + number of sideways neighboring blocks.
        /// This is due to the fact that the PeBi grid system implemented depends on flat blocks.
        /// So, each block will have a single top and a single bottom block except if it is a top/bottom boundary block.
        /// </remarks>
        public double[] deltaXList;

        /// <summary>
        /// An <see cref="Array"/> of the sides boundary lengths of the block.
        /// </summary>
        /// <remarks>
        /// The size the array is equal to 2 + number of sideways neighboring blocks.
        /// This is due to the fact that the PeBi grid system implemented depends on flat blocks.
        /// So, each block will have a single top and a single bottom block except if it is a top/bottom boundary block.
        /// </remarks>
        public double[] boundaryLengthList;

        /// <summary>
        /// The height of the block.
        /// </summary>
        public double height;

        /// <summary>
        /// The bulk volume of the block.
        /// </summary>
        public double bulkVolume;

        /// <summary>
        /// Pore volume array.
        /// </summary>
        /// <remarks>
        /// The size of this array depends on the number of steps memory used.
        /// </remarks>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double[] Vp;

        #endregion

        #region Well data
        /// <summary>
        /// The type of a block.
        /// </summary>
        /// <remarks>
        /// This will indicate if it's a normall block, or if it cintains a well.
        /// </remarks>
        /// <seealso cref="Global.BlockType"/>
        public Global.BlockType type;

        #endregion

        #region Miscellaneous

        /// <summary>
        /// list of the blocks transmissibilities, beginning with the top and bottom ones.
        /// </summary>
        /// <remarks>
        /// <para>This is the geometric transmissibility factor between two blocks. It does not change with time.</para>
        /// <para>The size of this array is equal to 2 + number of sideways neighboring blocks.</para>
        /// </remarks>
        /// <seealso cref="Transmissibility.Calculate(BaseBlock, BaseBlock, double)"/>
        public double[] transmissibility_list;


        /// <summary>
        /// These variables are used to save repeated calculations unnecessarily.
        /// </summary>
        /// <remarks>
        /// <para>To reduce clutter, these values are calculated once and stored for subsequent reuse.</para>
        /// <para>Unlike <see cref="transmissibility_list"/>, this is the final pressure and saturation dependent transmissibility.</para>
        /// </remarks>
        public double[] transmissibility_terms_oil, transmissibility_terms_water, transmissibility_terms_gas;

        /// <summary>
        /// These variables are used to save repeated calculations unnecessarily.
        /// </summary>
        public double accumulation_term_oil, accumulation_term_water, accumulation_term_gas;

        #endregion

        public double Depth;

        /// <summary>
        /// Initializes a new instance of the <see cref="BaseBlock"/> class.
        /// </summary>
        public BaseBlock()
        {
            int steps_memory = Global.STEPS_MEMORY;

            // data.porosity_calculator
            this.porosity = new double[steps_memory];

            // fluid
            this.Bo = new double[steps_memory]; this.Bg = new double[steps_memory]; this.Bw = new double[steps_memory];

            this.So = new double[steps_memory]; this.Sg = new double[steps_memory]; this.Sw = new double[steps_memory];

            this.Kro = new double[steps_memory]; this.Krg = new double[steps_memory]; this.Krw = new double[steps_memory];

            this.viscosityOil = new double[steps_memory]; this.viscosityGas = new double[steps_memory]; this.viscosityWater = new double[steps_memory];

            this.Rso = new double[steps_memory];

            this.P = new double[steps_memory];

            this.Pcow = new double[steps_memory];

            this.Pcgo = new double[steps_memory];

            // volumetric
            this.Vp = new double[steps_memory];

            // well
            this.type = Global.BlockType.NormalBlock;
        }


        public double GetPw(int P_time_level, int S_time_level)
        {
            return P[P_time_level] - Pcow[S_time_level];
        }

        public double GetPg(int P_time_level, int S_time_level)
        {
            return P[P_time_level] + Pcgo[S_time_level];
        }
    }
}