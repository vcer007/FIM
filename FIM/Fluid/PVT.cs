using FIM.Core;
using System;
using System.Linq;

namespace FIM.Fluid
{
    //Method Name: lookUp
    //Objectives: this is a private (internal method for the internal use of this class only) method that makes table lookup interpolation
    //Inputs: two variables representing the X and Y arrays "The two columns of the table" and the Y value
    //Outputs: a variable that represents the X value corresponding to the particular Y value

    //Method Name: extrapolate
    //Objectives: this is a private (internal method for the internal use of this class only) method that makes table lookup and extrapolation
    //Inputs: two variables representing the X and Y arrays "The two columns of the table" and the Y value
    //Outputs: a variable that represents the X value corresponding to the particular Y value

    //Method Name: getPVT "where PVT may be OilFVF, GasViscosity, WaterDensity, ..."
    //Objectives: calculate the corresponding PVT value at a certain pressure
    //Inputs: a variable representing the value of the pressure
    //Outputs: the corresponding PVT value
    class PVT
    {
        private double bubble_point_pressure;
        private double[][] oil_data, oil_us_data, water_data, water_us_data, gas_data;

        //Constructor Name: PVT
        //Objectives: Public constructor of the class. This is used only once to get the PVT data stored within the PVT class instance
        //Inputs: arrays representing the PVT values
        public PVT(double[][] oil_data = null, double[][] oil_us_data = null, double[][] water_data = null, double[][] water_us_data = null, double[][] gas_data = null, double bubble_point = 14.7)
        {
            this.oil_data = oil_data;
            this.oil_us_data = oil_us_data;

            this.water_data = water_data;
            this.water_us_data = water_us_data;

            this.gas_data = gas_data;

            this.bubble_point_pressure = bubble_point;
        }

        

        // The publicly accessible methods.
        public double getFVF(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return getWaterFVF(pressure);
                case Global.Phase.Oil:
                    return getOilFVF(pressure);
                case Global.Phase.Gas:
                    return getGasFVF(pressure);
                default:
                    return 1;
            }
        }

        public double getViscosity(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return getWaterViscosity(pressure);
                case Global.Phase.Oil:
                    return getOilViscosity(pressure);
                case Global.Phase.Gas:
                    return getGasViscosity(pressure);
                default:
                    return 1;
            }
        }

        public double getDensity(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return getWaterDensity(pressure);
                case Global.Phase.Oil:
                    return getOilDensity(pressure);
                case Global.Phase.Gas:
                    return getGasDensity(pressure);
                default:
                    return 1;
            }
        }

        public double getRs(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return getRsw(pressure);
                case Global.Phase.Oil:
                    return getRso(pressure) / Global.a;
                case Global.Phase.Gas:
                    return 1;
                default:
                    return 1;
            }
        }



        //Internal helper method to do the table lookups and interpolations
        private double lookUp(double[] data_y, double[] data_x, double y)
        {
            double y1, y2, x1, x2, x = -1;
            for (int i = 0; i < data_y.Length; i++)
            {
                //if an exact match exists, return the table value. No interpolation is required
                if (data_y[i] == y)
                {
                    return data_x[i];
                }

                //Look for the first element in the data_y array that is smaller than y
                if (data_y[i] > y)
                {
                    y1 = data_y[i - 1]; x1 = data_x[i - 1];
                    y2 = data_y[i]; x2 = data_x[i];
                    x = x1 + (y - y1) / ((y2 - y1) / (x2 - x1));
                    break;
                }
            }

            return x;
        }

        //Internal helper method for extrapolation
        private double extrapolate(double[] data_y, double[] data_x, double y)
        {
            double y1, y2, x1, x2, x = -1;

            y1 = data_y[data_y.Length - 2]; x1 = data_x[data_x.Length - 2];
            y2 = data_y.Last(); x2 = data_x.Last();

            x = x1 + (y - y1) / ((y2 - y1) / (x2 - x1));

            return x >= 0 ? x : 0;
        }


        // private methods used for the internal plumbing.
        private double getOilFVF(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                Y = oil_us_data[0]; X = oil_us_data[1];
                return extrapolate(Y, X, pressure);
            }
            else
            {
                Y = oil_data[0]; X = oil_data[1];
                return lookUp(Y, X, pressure);
            }

        }
        private double getWaterFVF(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                Y = water_us_data[0]; X = water_us_data[1];
                return extrapolate(Y, X, pressure);
            }
            else
            {
                Y = water_data[0]; X = water_data[1];
                return lookUp(Y, X, pressure);
            }
        }
        private double getGasFVF(double pressure)
        {
            double[] Y, X;
            Y = gas_data[0]; X = gas_data[1];

            double FVF = lookUp(Y, X, pressure);

            return FVF > 0 ? FVF : getGasFVF(gas_data[0].Last());
        }

        private double getOilViscosity(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                Y = oil_us_data[0]; X = oil_us_data[2];
                return extrapolate(Y, X, pressure);
            }
            else
            {
                Y = oil_data[0]; X = oil_data[2];
                return lookUp(Y, X, pressure);
            }
        }
        private double getWaterViscosity(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                Y = water_us_data[0]; X = water_us_data[2];
                return extrapolate(Y, X, pressure);
                //return lookUp(Y, X, pressure);

            }
            else
            {
                Y = water_data[0]; X = water_data[2];
                return lookUp(Y, X, pressure);
            }
        }
        private double getGasViscosity(double pressure)
        {
            double[] Y, X;
            Y = gas_data[0]; X = gas_data[2];

            double FVF = lookUp(Y, X, pressure);

            return pressure <= gas_data[0].Last() ? FVF : extrapolate(Y, X, pressure);
        }

        private double getOilDensity(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                Y = oil_us_data[0]; X = oil_us_data[3];
                return extrapolate(Y, X, pressure);
            }
            else
            {
                Y = oil_data[0]; X = oil_data[3];
                return lookUp(Y, X, pressure);
            }

        }
        private double getWaterDensity(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                Y = water_us_data[0]; X = water_us_data[3];
                return extrapolate(Y, X, pressure);
            }
            else
            {
                Y = water_data[0]; X = water_data[3];
                return lookUp(Y, X, pressure);
            }
        }
        private double getGasDensity(double pressure)
        {
            double[] Y, X;
            Y = gas_data[0]; X = gas_data[3];

            double Density = lookUp(Y, X, pressure);

            return pressure <= gas_data[0].Last() ? Density : extrapolate(Y, X, pressure);
        }

        private double getRso(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                return oil_us_data[4][0];
            }
            else
            {
                Y = oil_data[0]; X = oil_data[4];
                return lookUp(Y, X, pressure);
            }
        }
        private double getRsw(double pressure)
        {
            double[] Y, X;

            if (pressure >= bubble_point_pressure)
            {
                return water_us_data[4][0];
            }
            else
            {
                Y = water_data[0]; X = water_data[4];
                return lookUp(Y, X, pressure);
            }
        }

        //Method name: chord_slope_FVF
        //Objectives: calculates the FVF property at a future point pressure "P at time n+1" using the value of compressibility
        //Inputs: the value of the compressibility and the new pressure to calculate the pressure difference
        //Outputs: the new value of the FVF propery "FVF at time n+1"
        public static double chord_slope_FVF(double C, double new_pressure)
        {
            return 1 / (1 + C * (new_pressure - 14.7));
        }

        //Method Name: initializePressure
        //Objectives: initializing the pressure for problems considering the gravity effect
        //Inputs: there are two overloads "versions" of this method. one that takes an input for fluid compressibility "C" used for incompressible and slightly-compressible fluids
        //Inputs: the other form takes no input for compressibility, but will use the PVT values provided for density calculations
        //Outputs: an array of the values of the pressure for each depth
        public static double[] initializePressure(double[] depth_array, double pressure, double pressure_measurement_depth, double density, double C, double FVF)
        {
            double convergence_criterion = 0.001;
            int length = depth_array.Length;
            double[] pressure_array = new double[length];
            double[] old_pressure = new double[length];
            double[] density_array = new double[length];
            //
            //Initialize density-array
            for (int i = 0; i < length; i++)
            {
                depth_array[i] = density;
            }

            double density_standard = density * FVF;
            double new_FVF;

            //if the fluid is incompressible, then no-iterations are needed
            if (C == 0)
            {
                for (int i = 0; i < length; i++)
                {
                    pressure_array[i] = pressure + density * (depth_array[i] - pressure_measurement_depth);
                }
            }
            //iterate for the value of density for slightly-compressible fluids
            else
            {
                //First iteration
                for (int i = 0; i < length; i++)
                {
                    pressure_array[i] = pressure + density_array[i] * (depth_array[i] - pressure_measurement_depth);
                }
                //update density-array
                for (int i = 0; i < length; i++)
                {
                    new_FVF = chord_slope_FVF(C, pressure_array[i]);
                    density_array[i] = density_standard / new_FVF;
                }

                //From the second iteration on, check for convergence
                while (checkPressureConvergence(pressure_array, old_pressure) >= convergence_criterion)
                {
                    old_pressure = pressure_array.Clone() as double[];

                    for (int i = 0; i < length; i++)
                    {
                        pressure_array[i] = pressure + density_array[i] * (depth_array[i] - pressure_measurement_depth);
                    }
                    //update density-array
                    for (int i = 0; i < length; i++)
                    {
                        new_FVF = chord_slope_FVF(C, pressure_array[i]);
                        density_array[i] = density_standard / new_FVF;
                    }
                }
            }

            return pressure_array;
        }

        //Method Name: checkPressureConvergence
        //Objectives: used to check if the pressure solution converges through the non-linear iterations
        //Inputs: two arrays with the values of the old and new grid blocks pressures
        //Outputs: the maximum change in pressure that occurs in any of the grid blocks between two successive non-linear iterations
        //The output is compared against a specified criterion of convergence
        public static double checkPressureConvergence(double[] new_P, double[] old_P)
        {
            int length = new_P.Length;
            double[] temp = new double[length];
            for (int i = 0; i < length; i++)
            {
                temp[i] = Math.Abs(new_P[i] - old_P[i]) / old_P[i];
            }

            double max = temp.Max();
            return max;
        }

        // uses a Newton-Rapshon iteration to calculate the natural gas Z factor based on Dranchuk and Abu-Kassem equations
        public static double calculateZ(double Pc, double Tc, double P, double T, double z_initial)
        {
            double[] A = new double[] { 0.3265, -1.07, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.721 };

            double z, density_reduced, Pr, Tr;

            Pr = P / Pc;
            Tr = (T + 460) / Tc;

            z = z_initial;

            double C0 = (A[0] + A[1] / Tr + A[2] / Math.Pow(Tr, 3) + A[3] / Math.Pow(Tr, 4) + A[4] / Math.Pow(Tr, 5));
            double C1 = (A[5] + A[6] / Tr + A[7] / Math.Pow(Tr, 2));
            double C2 = -1 * A[8] * (A[6] / Tr + A[7] / Math.Pow(Tr, 2));

            density_reduced = 0.27 * Pr / (z * Tr);

            double f, f_dash;

            int counter = 0;
            do
            {
                f = 1
                    + C0 * density_reduced
                    + C1 * Math.Pow(density_reduced, 2)
                    + C2 * Math.Pow(density_reduced, 5)
                    + A[9] * (1 + A[10] * Math.Pow(density_reduced, 2)) * Math.Pow(density_reduced, 2) / Math.Pow(Tr, 3)
                    * Math.Exp(-1 * A[10] * Math.Pow(density_reduced, 2))
                    - 0.27 * Pr / (density_reduced * Tr)
                ;

                f_dash = C0
                    + 2 * C1 * density_reduced
                    + 5 * C2 * Math.Pow(density_reduced, 4)
                    + A[9] * (1 + A[10] * Math.Pow(density_reduced, 2)) * Math.Pow(density_reduced, 2) / Math.Pow(Tr, 3)
                    * Math.Exp(-1 * A[10] * Math.Pow(density_reduced, 2))
                    * -2 * A[10] * density_reduced
                    + Math.Exp(-1 * A[10] * Math.Pow(density_reduced, 2))
                    * A[9] * (2 + 4 * A[10] * Math.Pow(density_reduced, 2)) * density_reduced / Math.Pow(Tr, 3)
                    + 0.27 * Pr / (Tr * Math.Pow(density_reduced, 2))
                    ;

                density_reduced = density_reduced - f / f_dash;
                counter += 1;
            } while (counter < 5);

            return 0.27 * Pr / (density_reduced * Tr);
        }

        // calculate natural gas density 
        public static double calculateGasDensity(double pressure, double molecular_weight, double z_factor, double temperature)
        {
            return pressure * molecular_weight / (z_factor * (temperature + 460) * 10.73);
        }
    }
}
