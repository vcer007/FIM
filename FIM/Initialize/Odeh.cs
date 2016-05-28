using FIM.Core;
using FIM.Fluid;
using FIM.Misc;
using FIM.Rock;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Initialize
{
    class Odeh
    {
        public static SimulationData initiaize()
        {
            SimulationData simulation_data;
            PVT pvt;
            Kr kr;
            Porosity porosity;

            initializeFluidData(out pvt, out kr);
            intializePorosity(out porosity);

            initializeGrid(out simulation_data, pvt, kr);
            // initialize "n1" time level intitial guess
            initializen1(simulation_data);

            initializeTransmissibilities(simulation_data);

            simulation_data.phases = new Global.Phase[] { Global.Phase.Water, Global.Phase.Oil, Global.Phase.Gas };
            simulation_data.solubleGasPresent = true;

            simulation_data.time_stpe = 10;

            for (int i = 0; i < simulation_data.grid.Length; i++)
            {
                simulation_data.grid[i].updateProperties(pvt, kr, porosity, 4800, 0.12, 0.88, 0);
            }

            return simulation_data;
        }

        private static void initializen1(SimulationData simulation_data)
        {
            BaseBlock block;

            for (int i = 0; i < simulation_data.grid.Length; i++)
            {
                block = simulation_data.grid[i];

                // rock
                block.porosity[1] = block.porosity[0];

                // fluid
                block.Bo[1] = block.Bo[0]; block.Bg[1] = block.Bg[0]; block.Bw[1] = block.Bw[0];
                block.So[1] = block.So[0]; block.Sg[1] = block.Sg[0]; block.Sw[1] = block.Sw[0];
                block.Kro[1] = block.Kro[0]; block.Krg[1] = block.Krg[0]; block.Krw[1] = block.Krw[0];
                block.viscosity_oil[1] = block.viscosity_oil[0]; block.viscosity_gas[1] = block.viscosity_gas[0]; block.viscosity_water[1] = block.viscosity_water[0];
                block.Rso[1] = block.Rso[0];

                block.P[1] = block.P[0];  block.Po[1] = block.Po[0]; block.Pg[1] = block.Pg[0]; block.Pw[1] = block.Pw[0];

                // volumetric
                block.Vp[1] = block.Vp[0];

                // well
                block.flow_rate[1] = block.flow_rate[0]; block.BHP[1] = block.BHP[0]; block.J[1] = block.J[0];
                block.q_oil[1] = block.q_oil[0]; block.q_gas[1] = block.q_gas[0]; block.q_water[1] = block.q_water[0];

            }
        }

        private static void intializePorosity(out Porosity porosity)
        {
            double Cf = 0.000003;
            double porosity_ref = 0.3;
            double porosity_pressure_ref = 14.7;

            porosity = new Porosity(Cf, porosity_ref, porosity_pressure_ref);
        }

        private static void initializeGrid(out SimulationData simulation_data, PVT pvt, Kr kr)
        {
            int x = 10, y = 10, z = 3;

            double porosity = 0.3;
            double[][] permeability = new double[3][];
            permeability[0] = new double[] { 50, 50, 500};
            permeability[1] = new double[] { 50, 25, 50 };
            permeability[2] = new double[] { 25, 25, 200 };

            double Sw = 0.12, So = 0.88, Sg = 0;
            double pressure = 4800;

            double delta_x = 1000, delta_y = 1000;
            double[] h = new double[] { 20, 30, 50 };

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
                        grid[counter].neighbour_blocks_indices = temp.ToArray();
                        grid[counter].index = Rectangular.xyzToNatural(x, y, z, i, j, k);
                        grid[counter].layer = k;

                        //
                        grid[counter].transmissibility_terms_oil = new double[temp.Count];
                        grid[counter].transmissibility_terms_water = new double[temp.Count];
                        grid[counter].transmissibility_terms_gas = new double[temp.Count];

                        temp.Clear();
                        counter += 1;
                    }
                }
            }


            // initialize rock properties
            for (int i = 0; i < grid.Length; i++)
            {
                size = grid[i].neighbour_blocks_indices.Length;

                grid[i].porosity[0] = porosity;
                grid[i].permeability_list = new double[size];

                grid[i].permeability_list[0] = permeability[grid[i].layer][0];
                grid[i].permeability_list[1] = permeability[grid[i].layer][1];

                for (int a = 2; a < size; a++)
                {
                    grid[i].permeability_list[a] = permeability[grid[i].layer][2];
                }
            }

            // intialize fluid properties
            for (int i = 0; i < grid.Length; i++)
            {
                grid[i].Sw[0] = Sw; grid[i].So[0] = So; grid[i].Sg[0] = Sg;

                grid[i].Kro[0] = kr.getKr(Global.Phase.Oil, Sg);
                grid[i].Krw[0] = kr.getKr(Global.Phase.Water, Sg);
                grid[i].Krg[0] = kr.getKr(Global.Phase.Gas, Sg);

                grid[i].Bo[0] = pvt.getFVF(Global.Phase.Oil, pressure);
                grid[i].Bw[0] = pvt.getFVF(Global.Phase.Water, pressure);
                grid[i].Bg[0] = pvt.getFVF(Global.Phase.Gas, pressure);

                grid[i].viscosity_oil[0] = pvt.getViscosity(Global.Phase.Oil, pressure);
                grid[i].viscosity_water[0] = pvt.getViscosity(Global.Phase.Water, pressure);
                grid[i].viscosity_gas[0] = pvt.getViscosity(Global.Phase.Gas, pressure);

                grid[i].Rso[0] = pvt.getRs(Global.Phase.Oil, pressure);

                grid[i].P[0] = pressure;
            }

            // intialize volumetric data
            for (int i = 0; i < grid.Length; i++)
            {
                double height = h[grid[i].layer];
                grid[i].h = height;
                grid[i].delta_x_list = new double[6];

                grid[i].delta_x_list[0] = height;
                grid[i].delta_x_list[1] = height;

                for (int a = 2; a < grid[i].delta_x_list.Length; a++)
                {
                    grid[i].delta_x_list[a] = delta_x;
                }

                // areas
                grid[i].area_list = new double[6];

                grid[i].area_list[0] = delta_x * delta_y;
                grid[i].area_list[1] = delta_x * delta_y;

                for (int a = 2; a < grid[i].area_list.Length; a++)
                {
                    grid[i].area_list[a] = grid[i].h * delta_x;
                }

                grid[i].bulk_volume = grid[i].h * delta_x;

                grid[i].Vp[0] = grid[i].bulk_volume * porosity;
            }

            simulation_data = new SimulationData(x, y, z, grid);
        }

        private static void initializeFluidData(out PVT pvt, out Kr kr)
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
            kr = new Kr(Kr_data);
        }

        private static void initializeTransmissibilities(SimulationData simulation_data)
        {
            for (int i = 0; i < simulation_data.grid.Length; i++)
            {
                simulation_data.grid[i].transmissibility_list = new double[simulation_data.grid[i].neighbour_blocks_indices.Length];

                for (int a = 0; a < simulation_data.grid[i].neighbour_blocks_indices.Length; a++)
                {
                    if (simulation_data.grid[i].neighbour_blocks_indices[a] != -1)
                    {
                        simulation_data.grid[i].transmissibility_list[a] = Transmissibility.calculate(simulation_data.grid[i], simulation_data.grid[simulation_data.grid[i].neighbour_blocks_indices[a]]);
                    }
                    else
                    {
                        simulation_data.grid[i].transmissibility_list[a] = 0;
                    }
                }
            }
        }
    }
}
