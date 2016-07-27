using FIM.Core;
using FIM.FluidData;
using FIM.Misc;
using FIM.RockData;
using FIM.Well;
using FIM.Extensions;

using System.Collections.Generic;

/// <summary>
/// This name space contains the classes and methods used for intializing models..
/// </summary>
namespace FIM.Initialization
{
    static class Model
    {
        public static SimulationData initiaize()
        {
            SimulationData simulation_data;
            PVT pvt;
            SCAL kr;
            PorosityCalculator porosity;

            initializeFluidData(out pvt, out kr);
            intializePorosityCalculator(out porosity);

            initializeGrid(out simulation_data, pvt, kr);

            initializeTransmissibilities(simulation_data);

            simulation_data.phases = new Global.Phase[] { Global.Phase.Oil, Global.Phase.Gas, Global.Phase.Water };
            simulation_data.solubleGasPresent = true;
            simulation_data.solutionProcedure = Global.SolutionProcedure.FullyImplicit;

            simulation_data.originalTimeStep = 1;
            simulation_data.timeStepSlashingFactor = 0.5;
            simulation_data.endingTime = 10 * 365;

            simulation_data.originalRelaxationFactor = 1;
            simulation_data.minimumRelaxation = 0.5;
            simulation_data.relaxationFactorIncrement = -0.1;
            simulation_data.relaxationFactor = simulation_data.originalRelaxationFactor;

            simulation_data.maximumNonLinearIterations = 25;

            simulation_data.maximumConvergenceErrorRatio = 0.5;

            simulation_data.pvt = pvt;
            simulation_data.scal = kr;
            simulation_data.porosityCalculator = porosity;

            initializeWells(simulation_data);


            for (int i = 0; i < simulation_data.grid.Length; i++)
            {
                simulation_data.grid[i].UpdateProperties(simulation_data, 4800, 0.12, 0, 0.88, 0);
            }

            simulation_data.MBE_Tolerance = 1;


            return simulation_data;
        }

        private static void initializeFluidData(out PVT pvt, out SCAL kr)
        {
            double[][] oil, oil_us, water, water_us, gas;

            oil = new double[5][]; oil_us = new double[5][]; water = new double[5][]; water_us = new double[5][]; gas = new double[4][];

            oil[0] = new double[] { 14.7, 264.7, 514.7, 1014.7, 2014.7, 2514.7, 3014.7, 4014.7, 5014.7, 9014.7 };
            oil[1] = new double[] { 1.062, 1.15, 1.207, 1.295, 1.435, 1.5, 1.565, 1.695, 1.827, 2.357 };
            oil[2] = new double[] { 1.04, 0.975, 0.91, 0.83, 0.695, 0.641, 0.594, 0.51, 0.449, 0.203 };
            oil[3] = new double[] { 46.244, 43.544, 42.287, 41.004, 38.995, 38.304, 37.781, 37.046, 36.424, 34.482 };
            oil[4] = new double[] { 1, 90.5, 180, 371, 636, 775, 930, 1270, 1618, 2984 };

            oil_us[0] = new double[] { 4014.7, 9014.7 };
            oil_us[1] = new double[] { 1.695, 1.579 };
            oil_us[2] = new double[] { 0.51, 0.74 };
            oil_us[3] = new double[] { 37.046, 39.768 };
            oil_us[4] = new double[] { 1270, 1270 };

            water[0] = new double[] { 14.7, 264.7, 514.7, 1014.7, 2014.7, 2514.7, 3014.7, 4014.7, 5014.7, 9014.7 };
            water[1] = new double[] { 1.041, 1.0403, 1.0395, 1.038, 1.035, 1.0335, 1.032, 1.029, 1.0258, 1.013 };
            water[2] = new double[] { 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31 };
            water[3] = new double[] { 62.238, 62.283, 62.328, 62.418, 62.599, 62.69, 62.781, 62.964, 63.16, 63.959 };
            water[4] = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            water_us[0] = new double[] { 4014.7, 9014.7 };
            water_us[1] = new double[] { 1.029, 1.013 };
            water_us[2] = new double[] { 0.31, 0.31 };
            water_us[3] = new double[] { 62.964, 63.959 };
            water_us[4] = new double[] { 0, 0 };

            gas[0] = new double[] { 14.7, 264.7, 514.7, 1014.7, 2014.7, 2514.7, 3014.7, 4014.7, 5014.7, 9014.7 };
            gas[1] = new double[] { 0.166666, 0.012093, 0.006274, 0.003197, 0.001614, 0.001294, 0.00108, 0.000811, 0.000649, 0.000386 };
            gas[2] = new double[] { 0.008, 0.0096, 0.0112, 0.014, 0.0189, 0.0208, 0.0228, 0.0268, 0.0309, 0.047 };
            gas[3] = new double[] { 0.0647, 0.8916, 1.7185, 3.3727, 6.8806, 8.3326, 9.9837, 13.2952, 16.6139, 27.948 };

            double[][] Kr_data = new double[4][];
            Kr_data[0] = new double[] { 0, 0.001, 0.02, 0.05, 0.12, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.85, 1 };
            Kr_data[1] = new double[] { 0, 0, 0, 0.005, 0.025, 0.075, 0.125, 0.19, 0.41, 0.6, 0.72, 0.87, 0.94, 0.98, 1 };
            Kr_data[2] = new double[] { 1, 1, 0.997, 0.98, 0.7, 0.35, 0.2, 0.09, 0.02, 0.01, 0.001, 0.0001, 0, 0, 0 };
            Kr_data[3] = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            pvt = new PVT(oil, oil_us, water, water_us, gas, 4014.7);
            kr = new SCAL(Kr_data);
        }

        private static void intializePorosityCalculator(out PorosityCalculator porosity)
        {
            double Cf = 0.000003;
            double porosity_ref = 0.3;
            double porosity_pressure_ref = 14.7;

            porosity = new PorosityCalculator(Cf, porosity_ref, porosity_pressure_ref);
        }

        private static void initializeGrid(out SimulationData simulation_data, PVT pvt, SCAL kr)
        {
            int x = 2, y = 1, z = 1;

            double porosity = 0.3;
            double[][] permeability = new double[3][];
            permeability[0] = new double[] { 50, 50, 500 };
            permeability[1] = new double[] { 50, 25, 50 };
            permeability[2] = new double[] { 25, 25, 200 };

            double delta_x = 1000, delta_y = 1000;
            double[] h = new double[] { 20, 30, 50};

            int size = x * y * z;

            BaseBlock[] grid = new BaseBlock[size];

            List<int> temp = new List<int>();

            // intialize rectangular blocks neighbours list
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
                        grid[counter].neighbourBlocksIndices = temp.ToArray();
                        grid[counter].index = Rectangular.xyzToNatural(x, y, z, i, j, k);
                        grid[counter].layer = k;

                        grid[counter].boundaryLengthList = new double[temp.Count];
                        for (int a = 0; a < temp.Count; a++)
                        {
                            // for the Odeh problem, delta_x = delta_y = 1000 for all blocks.
                            grid[counter].boundaryLengthList[a] = delta_x;
                        }

                        //
                        grid[counter].transmissibility_terms_oil = new double[temp.Count];
                        grid[counter].transmissibility_terms_water = new double[temp.Count];
                        grid[counter].transmissibility_terms_gas = new double[temp.Count];

                        temp.Clear();
                        counter += 1;
                    }
                }
            }


            // initialize porosity_calculator properties
            for (int i = 0; i < grid.Length; i++)
            {
                size = grid[i].neighbourBlocksIndices.Length;

                grid[i].porosity[0] = porosity;
                grid[i].permeability = new double[size];

                grid[i].permeability[0] = permeability[grid[i].layer][0];
                grid[i].permeability[1] = permeability[grid[i].layer][1];

                for (int a = 2; a < size; a++)
                {
                    grid[i].permeability[a] = permeability[grid[i].layer][2];
                }
            }

            // intialize volumetric data
            for (int i = 0; i < grid.Length; i++)
            {
                double height = h[grid[i].layer];
                grid[i].height = height;
                grid[i].deltaXList = new double[6];

                grid[i].deltaXList[0] = height;
                grid[i].deltaXList[1] = height;

                for (int a = 2; a < grid[i].deltaXList.Length; a++)
                {
                    grid[i].deltaXList[a] = delta_x;
                }

                // areas
                grid[i].areaList = new double[6];

                grid[i].areaList[0] = delta_x * delta_y;
                grid[i].areaList[1] = delta_x * delta_y;

                for (int a = 2; a < grid[i].areaList.Length; a++)
                {
                    grid[i].areaList[a] = grid[i].height * delta_x;
                }

                grid[i].bulkVolume = grid[i].height * delta_x * delta_x;

                grid[i].Vp[0] = grid[i].bulkVolume * porosity;
            }

            simulation_data = new SimulationData(grid);
        }

        private static void initializeTransmissibilities(SimulationData simulation_data)
        {
            for (int i = 0; i < simulation_data.grid.Length; i++)
            {
                simulation_data.grid[i].transmissibility_list = new double[simulation_data.grid[i].neighbourBlocksIndices.Length];

                for (int a = 0; a < simulation_data.grid[i].neighbourBlocksIndices.Length; a++)
                {
                    if (simulation_data.grid[i].neighbourBlocksIndices[a] != -1)
                    {
                        simulation_data.grid[i].transmissibility_list[a] = Transmissibility.calculate(simulation_data.grid[i], simulation_data.grid[simulation_data.grid[i].neighbourBlocksIndices[a]]);
                    }
                    else
                    {
                        simulation_data.grid[i].transmissibility_list[a] = 0;
                    }
                }
            }
        }

        private static void initializeWells(SimulationData data)
        {
            int[] well_indices = new int[] { 0};

            // well data
            BaseWell[] wells = new BaseWell[1];

            wells[0] = new BaseWell(data, well_indices[0], Global.WellType.Production, Global.WellControl.OilRate, Global.WellRateCalculation.Explicit, 0.25, 0, 500, 1000);

            data.wells = wells;

        }

    }
}
