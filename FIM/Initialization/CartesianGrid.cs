using FIM.Core;
using FIM.Misc;

using System.Collections.Generic;
using System;
using System.Linq;

namespace FIM.Initialization
{
    // this class is solely used for a grid structure that is similar to the SPE1 comparative problem model.
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

        public BaseBlock[] Initialize()
        {
            // initialize the grid.
            // assign neighbor blocks indices and blocks arrays sizes.
            BaseBlock[] grid = InitializeGrid();

            // assign permeabilities
            // note that porosities "at different time levels" are calculated when updating blocks properties automaticaly.
            // this is due to the fact that porosity is pressure dependent. So it's calculated by knowing the rock's compressibility.
            InitializeGridPermeabilities(grid);

            InitializeVolumetrics(grid);

            InitializeTransmissibilities(grid);

            return grid;
        }

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

        private void InitializeGridPermeabilities(BaseBlock[] grid)
        {
            for (int i = 0; i < grid.Length; i++)
            {
                // assign each grid block's permeabilities.
                grid[i].permeability = Enumerable.Repeat(layersPermeabilities[grid[i].layer], grid[i].neighborBlocksIndices.Length).ToArray();
            }
        }

        public BaseBlock[] InitializeGrid()
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

                        InitializeBlockArrays(grid[counter], temp.Count);

                        InitializeBlockDimensions(i, j, k, grid[counter]);

                        temp.Clear();
                        counter += 1;
                    }
                }
            }

            return grid;
        }

        private void InitializeBlockArrays(BaseBlock block, int arraySize)
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

        private void InitializeBlockDimensions(int x, int y, int z, BaseBlock block)
        {
            block.height = layersHeights[z];

            // up
            block.deltaXList[0] = 0.5 * block.height;

            // down
            block.deltaXList[1] = block.deltaXList[0];

            // right
            block.deltaXList[2] = 0.5 * x_dimensionLengths[x];
            block.boundaryLengthList[0] = y_dimensionLengths[y];

            // left
            block.deltaXList[3] = block.deltaXList[2];
            block.boundaryLengthList[1] = block.boundaryLengthList[0];

            // north
            block.deltaXList[4] = 0.5 * y_dimensionLengths[y];
            block.boundaryLengthList[2] = x_dimensionLengths[x];

            // south
            block.deltaXList[5] = block.deltaXList[4];
            block.boundaryLengthList[3] = block.boundaryLengthList[2];

        }

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
