using FIM.Core;
using FIM.Misc;

using System.Collections.Generic;
using System.Linq;
using System;

namespace FIM.Initialization
{
    // this class is solely used for a grid structure that is similar to the SPE1 comparative problem model.
    // that is a cartesian grid.
    // although the simulator is perfectly capable to handle PeBi grids, cartesian grids are considered as a special case of a general PeBi grid.
    class CartesianGrid
    {
        // the grid x, y and z dimensions.
        int x, y, z;

        // permeabilities at each layer "z", at each direction "x, y, z".
        // a transmissibility multiplier is used to alter up/down directions transmissibilities.
        // it's assumed that an isotropic rock is present. Kx = Ky = Kz.
        double[] layersPermeabilities;

        public double[] transmissibilityMultipliers;

        // the heights of the layers "z".
        // the size of this array is equal to the number of blocks at the z direction "z".
        double[] layersHeights;

        // blocks longitudes "x direction block length" 
        // the size of this array is equal to the number of blocks at the x direction "x".
        double[] x_dimensionLengths;

        // blocks latitudes "y direction block length" 
        // the size of this array is equal to the number of blocks at the y direction "y".
        double[] y_dimensionLengths;

        /// <summary>
        /// Initializes a new instance of the <see cref="CartesianGrid"/> class.
        /// </summary>
        /// <remarks>
        /// This constructor is used when the blocks width and length "x and y dimension length" are different for each block along their respective direction.
        /// They are defines as arrays of length equal to the blocks number along each direction.
        /// </remarks>
        /// <param name="x">The x grid dimension.</param>
        /// <param name="y">The y grid dimension.</param>
        /// <param name="z">The z grid dimension.</param>
        /// <param name="layersPermeabilities">
        /// The layers permeabilities.
        /// Each layer is assumed to be isotropic; it has the same permeability regardless of direction.
        /// Multipliers <see cref="transmissibilityMultipliers"/> are used for modifying inter-layers transmissibilities.
        /// </param>
        /// <param name="x_dimensionLengths">The x_dimension lengths.</param>
        /// <param name="y_dimensionLengths">The y_dimension lengths.</param>
        /// <param name="layersHeights">The layers heights.</param>
        public CartesianGrid(int x, int y, int z, double[] layersPermeabilities, double[] x_dimensionLengths, double[] y_dimensionLengths, double[] layersHeights)
        {
            this.x = x;
            this.y = y;
            this.z = z;

            this.layersPermeabilities = layersPermeabilities;
            this.layersHeights = layersHeights;
            this.x_dimensionLengths = x_dimensionLengths;
            this.y_dimensionLengths = y_dimensionLengths;

        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CartesianGrid"/> class.
        /// </summary>
        /// <remarks>
        /// This constructor is used when the blocks width and length "x and y dimension length" are each a fixed single value.
        /// </remarks>
        /// <param name="x">The x grid dimension.</param>
        /// <param name="y">The y grid dimension.</param>
        /// <param name="z">The z grid dimension.</param>
        /// <param name="layersPermeabilities">The layers permeabilities.</param>
        /// <param name="x_dimensionLengths">The x_dimension length.</param>
        /// <param name="y_dimensionLengths">The y_dimension length.</param>
        /// <param name="layersHeights">The layers heights.</param>
        public CartesianGrid(int x, int y, int z, double[] layersPermeabilities, double x_dimensionLengths, double y_dimensionLengths, double[] layersHeights)
        {
            this.x = x;
            this.y = y;
            this.z = z;

            this.layersPermeabilities = layersPermeabilities;
            this.layersHeights = layersHeights;
            this.x_dimensionLengths = Enumerable.Repeat(x_dimensionLengths, x).ToArray();
            this.y_dimensionLengths = Enumerable.Repeat(y_dimensionLengths, y).ToArray();

        }


        /// <summary>
        /// Constructs the <see cref="BaseBlock"/> array and assign volumetric properties to the blocks.
        /// </summary>
        /// <returns> The <see cref="BaseBlock"/> array used as the grid in the <see cref="SimulationData"/></returns>
        /// <seealso cref="SimulationData.grid"/>
        public BaseBlock[] Initialize()
        {
            // initialize the grid.
            // assign neighbor blocks indices and blocks arrays sizes.
            BaseBlock[] grid = InitializeNeighbors();

            InitializeGridPermeabilities(grid);

            InitializeVolumetrics(grid);

            InitializeTransmissibilities(grid);

            InitializeDepths(grid);

            return grid;
        }

        private void InitializeDepths(BaseBlock[] grid)
        {
            for (int i = 0; i < grid.Length; i++)
            {
                grid[i].Depth = 0.5 * layersHeights[grid[i].layer] + layersHeights.Take(grid[i].layer).Sum();
            }
        }

        // an internal method to initialize blocks neghboring relationshps.
        private BaseBlock[] InitializeNeighbors()
        {
            int size = x * y * z;

            BaseBlock[] grid = new BaseBlock[size];

            List<int> temp = new List<int>();

            int counter = 0;
            for (int k = 0; k < z; k++)
            {
                for (int j = 0; j < y; j++)
                {
                    for (int i = 0; i < x; i++)
                    {
                        //up
                        if (k > 0) { temp.Add(Rectangular.xyzToNatural(x, y, z, i, j, k - 1)); } else { temp.Add(-1); };
                        //bottom
                        if (k < z - 1) { temp.Add(Rectangular.xyzToNatural(x, y, z, i, j, k + 1)); } else { temp.Add(-1); };
                        //right
                        if (i < x - 1) { temp.Add(Rectangular.xyzToNatural(x, y, z, i + 1, j, k)); } else { temp.Add(-1); };
                        //left
                        if (i > 0) { temp.Add(Rectangular.xyzToNatural(x, y, z, i - 1, j, k)); } else { temp.Add(-1); };
                        //north
                        if (j < y - 1) { temp.Add(Rectangular.xyzToNatural(x, y, z, i, j + 1, k)); } else { temp.Add(-1); };
                        //south
                        if (j > 0) { temp.Add(Rectangular.xyzToNatural(x, y, z, i, j - 1, k)); } else { temp.Add(-1); };

                        grid[counter] = new BaseBlock();
                        grid[counter].neighborBlocksIndices = temp.ToArray();
                        grid[counter].index = Rectangular.xyzToNatural(x, y, z, i, j, k);
                        grid[counter].layer = k;

                        AllocateBlockArrays(grid[counter], temp.Count);

                        InitializeBlockDimensions(i, j, k, grid[counter]);

                        temp.Clear();
                        counter += 1;
                    }
                }
            }

            return grid;
        }

        // allocates memory for each array in the block according to the number of neighbor blocks.
        // arraySize = the number of neighbor blocks.
        private void AllocateBlockArrays(BaseBlock block, int arraySize)
        {
            // boundary length list is only used for well equivalent radius calculations.
            // only sideways boundaries are considered.
            block.boundaryLengthList = new double[arraySize - 2];
            block.transmissibility_terms_oil = new double[arraySize];
            block.transmissibility_terms_water = new double[arraySize];
            block.transmissibility_terms_gas = new double[arraySize];
            block.deltaXList = new double[arraySize];
            block.areaList = new double[arraySize];
        }

        // initializes the blocks delta Xs.
        // i, j and k are the block's x, y and z location in the cartesian grid.
        private void InitializeBlockDimensions(int i, int j, int k, BaseBlock block)
        {
            block.height = layersHeights[k];

            // up
            block.deltaXList[0] = 0.5 * block.height;

            // down
            block.deltaXList[1] = block.deltaXList[0];

            // right
            block.deltaXList[2] = 0.5 * x_dimensionLengths[i];
            block.boundaryLengthList[0] = y_dimensionLengths[j];

            // left
            block.deltaXList[3] = block.deltaXList[2];
            block.boundaryLengthList[1] = block.boundaryLengthList[0];

            // north
            block.deltaXList[4] = 0.5 * y_dimensionLengths[j];
            block.boundaryLengthList[2] = x_dimensionLengths[i];

            // south
            block.deltaXList[5] = block.deltaXList[4];
            block.boundaryLengthList[3] = block.boundaryLengthList[2];

        }


        // assign permeabilities
        // note that porosities "at different time levels" are calculated when updating blocks properties automaticaly.
        // this is due to the fact that porosity is pressure dependent. So it's calculated by knowing the rock's compressibility.
        private void InitializeGridPermeabilities(BaseBlock[] grid)
        {
            for (int i = 0; i < grid.Length; i++)
            {
                // assign each grid block's permeabilities.
                grid[i].permeability = Enumerable.Repeat(layersPermeabilities[grid[i].layer], grid[i].neighborBlocksIndices.Length).ToArray();
            }
        }

        // assign boundary lengths, delta Xs, areas and block volumes.
        private void InitializeVolumetrics(BaseBlock[] grid)
        {
            for (int i = 0; i < grid.Length; i++)
            {
                // areas

                // top
                grid[i].areaList[0] = 4 * grid[i].deltaXList[2] * grid[i].deltaXList[4];

                // bottom
                grid[i].areaList[1] = grid[i].areaList[0];

                for (int a = 0; a < grid[i].boundaryLengthList.Length; a++)
                {
                    grid[i].areaList[a + 2] = grid[i].height * grid[i].boundaryLengthList[a];
                }

                grid[i].bulkVolume = grid[i].height * grid[i].areaList[0];
            }
        }

        // intializes the blocks geometrical transmissibility factors.
        private void InitializeTransmissibilities(BaseBlock[] grid)
        {
            for (int i = 0; i < grid.Length; i++)
            {
                grid[i].transmissibility_list = new double[grid[i].neighborBlocksIndices.Length];

                for (int a = 0; a < grid[i].neighborBlocksIndices.Length; a++)
                {
                    if (grid[i].neighborBlocksIndices[a] != -1)
                    {
                        grid[i].transmissibility_list[a] = Transmissibility.Calculate(grid[i], grid[grid[i].neighborBlocksIndices[a]]);
                    }
                    else
                    {
                        grid[i].transmissibility_list[a] = 0;
                    }


                    // apply transmissibility multipliers

                    // down
                    if (a == 1)
                    {
                        int layer = grid[i].layer;

                        grid[i].transmissibility_list[1] *= transmissibilityMultipliers[layer];

                    }
                    // up
                    else if (a == 0)
                    {
                        int layer = grid[i].layer - 1;

                        grid[i].transmissibility_list[0] *= layer == -1 ? 1 : transmissibilityMultipliers[layer];
                    }
                }
            }
        }
    }
}
