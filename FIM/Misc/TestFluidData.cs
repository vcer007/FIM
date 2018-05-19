using FIM.FluidData;

using System;

namespace FIM.Misc
{
    /// <summary>
    /// This class i made at an early stage of devloping the simulator to test the data accuracy and interpolation and extrapolation methods used.
    /// </summary>
    class TestFluidData
    {
        public static void test()
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

            PVT pvt = new PVT(oil, oil_us, water, water_us, gas, 4014.7);
            SCAL kr = new SCAL(Kr_data);

            for (double i = 0; i <= 1; i += 0.01)
            {
                //Console.WriteLine(i + ", " + kr.GetKr(Core.Global.Phase.Gas, i));
            }

            Console.ReadKey();
        }
    }
}
