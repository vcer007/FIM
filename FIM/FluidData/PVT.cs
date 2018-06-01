using FIM.Core;
using System;
using System.Linq;

namespace FIM.FluidData
{
    /// <summary>
    /// The Pressure Volume Temperature data storage class.
    /// </summary>
    /// <remarks>
    /// <para>This class acts as a container to store all input PVT data.</para>
    /// <para>It also contains methods for interpolation to get PVT data at different values of pressure.</para>
    /// </remarks>
    /// <seealso cref="SCAL"/>
    /// <seealso cref="RockData.PorosityCalculator"/>
    public class PVT
    {
        /// <summary>
        /// The bubble point pressure.
        /// </summary>
        public double bubblePointPressure, connateWaterSaturation;

        private double[] surfaceDensities, waterData;

        // these arrays are used to store the input PVT data.
        private double[][] oilData, oilUnderSaturatedData, gasData;



        /// <summary>
        /// Initializes a new instance of the <see cref="PVT"/> class.
        /// </summary>
        /// <remarks>
        /// This is used only once to get the PVT data stored within the PVT class instance
        /// </remarks>
        /// <param name="oilData">The saturated oil data.</param>
        /// <param name="oilUnderSaturatedData">The undersaturated oil data.</param>
        /// <param name="waterData">The water data.</param>
        /// <param name="waterUnderSaturatedData">The undersaturated water data.</param>
        /// <param name="gasData">The gas data.</param>
        /// <param name="bubblePointPressure">The bubble point pressure.</param>
        public PVT(double[][] oilData, double[][] oilUnderSaturatedData, double[] waterData, double[][] waterUnderSaturatedData, double[][] gasData, double bubblePointPressure)
        {
            this.oilData = oilData;
            this.oilUnderSaturatedData = oilUnderSaturatedData;

            this.waterData = waterData;

            this.gasData = gasData;

            this.bubblePointPressure = bubblePointPressure;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="PVT"/> class.
        /// </summary>
        /// <remarks>
        /// This is an empty constructor of a new PVT instance.
        /// The <see cref="Initialize(double[][], double[][], double[][], double[][], double[][], double)"/> method must be used to assign the data.
        /// </remarks>
        public PVT()
        {

        }

        internal double GetRv(double p)
        {
            var Y = oilData[0]; var X = oilData[3];
            var temp = LookUp(Y, X, p);
            return temp * Global.a / 1000; // cubic feet / scf
        }

        /// <summary>
        /// Initializes the member variables in a PVT instance.
        /// </summary>
        /// <remarks>
        /// This method is only used when a <see cref="PVT"/> instance is created using the empty constructor.
        /// The empty constructor is used to create a <see cref="PVT"/> object "instance" without specifying its data.
        /// This should be taken with so much care.
        /// </remarks>
        /// <param name="oilData">The oil data.</param>
        /// <param name="oilUnderSaturatedData">The oil under saturated data.</param>
        /// <param name="waterData">The water data.</param>
        /// <param name="waterUnderSaturatedData">The water under saturated data.</param>
        /// <param name="gasData">The gas data.</param>
        /// <param name="bubblePointPressure">The bubble point pressure.</param>
        public void Initialize(double[][] oilData = null, double[][] oilUnderSaturatedData = null, double[] waterData = null, double[][] gasData = null, double bubblePointPressure = 14.7, double[] surfaceDensities = null, double Swc = 0.12)
        {
            this.oilData = oilData;
            this.oilUnderSaturatedData = oilUnderSaturatedData;

            this.waterData = waterData;

            this.gasData = gasData;

            this.bubblePointPressure = bubblePointPressure;

            this.surfaceDensities = surfaceDensities;

            this.connateWaterSaturation = Swc;
        }

        // The publicly accessible methods.

        /// <summary>
        /// Gets the FVF of a certain phase.
        /// </summary>
        /// <param name="phase">The phase.</param>
        /// <param name="pressure">The pressure.</param>
        /// <returns>The value of FVF</returns>
        /// <seealso cref="Global.Phase"/>
        public double GetFVF(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return GetWaterFVF(pressure); // bbl / stb
                case Global.Phase.Oil:
                    return GetOilFVF(pressure); // bbl / stb
                case Global.Phase.Gas:
                    return GetGasFVF(pressure) * Global.a / 1000; // cubic feet per scf
                default:
                    return 1;
            }
        }

        /// <summary>
        /// Gets the viscosity of a certain phase.
        /// </summary>
        /// <param name="phase">The phase.</param>
        /// <param name="pressure">The pressure.</param>
        /// <returns>The value of viscosity</returns>
        /// <seealso cref="Global.Phase"/>
        public double GetViscosity(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return GetWaterViscosity(pressure);
                case Global.Phase.Oil:
                    return GetOilViscosity(pressure);
                case Global.Phase.Gas:
                    return GetGasViscosity(pressure);
                default:
                    return 1;
            }
        }

        /// <summary>
        /// Gets the density of a certain phase.
        /// </summary>
        /// <param name="phase">The phase.</param>
        /// <param name="pressure">The pressure.</param>
        /// <returns>The value of density</returns>
        /// <seealso cref="Global.Phase"/>
        //public double GetDensity(Global.Phase phase, double pressure)
        //{
        //    switch (phase)
        //    {
        //        case Global.Phase.Water:
        //            return GetWaterDensity(pressure);
        //        case Global.Phase.Oil:
        //            return GetOilDensity(pressure);
        //        case Global.Phase.Gas:
        //            return GetGasDensity(pressure);
        //        default:
        //            return 1;
        //    }
        //}

        /// <summary>
        /// Gets the solution Gas/Oil ratio.
        /// </summary>
        /// <param name="phase">The phase.</param>
        /// <param name="pressure">The pressure.</param>
        /// <returns>The value of RsO.</returns>
        public double GetRs(Global.Phase phase, double pressure)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return GetRsw(pressure);
                case Global.Phase.Oil:
                    return GetRso(pressure) / Global.a * 1000; // stb / bbl
                case Global.Phase.Gas:
                    return 1;
                default:
                    return 1;
            }
        }


        //Internal helper method to do the table lookups and interpolations
        //Method Name: lookUp
        //Objectives: this is a private (internal method for the internal use of this class only) method that makes table lookup interpolation
        //Inputs: two variables representing the X and Y arrays "The two columns of the table" and the Y value
        //Outputs: a variable that represents the X value corresponding to the particular Y value
        private double LookUp(double[] data_y, double[] data_x, double y)
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
        //Method Name: extrapolate
        //Objectives: this is a private (internal method for the internal use of this class only) method that makes table lookup and extrapolation
        //Inputs: two variables representing the X and Y arrays "The two columns of the table" and the Y value
        //Outputs: a variable that represents the X value corresponding to the particular Y value
        private double Extrapolate(double[] data_y, double[] data_x, double y)
        {
            double y1, y2, x1, x2, x = -1;

            int size = data_y.Length;

            y1 = data_y[data_y.Length - 2]; x1 = data_x[data_x.Length - 2];
            y2 = data_y[size - 1]; x2 = data_x[size - 1];

            x = x1 + (y - y1) / ((y2 - y1) / (x2 - x1));

            return x >= 0 ? x : 0;
        }


        // private methods used for the internal plumbing.

        private double GetOilFVF(double pressure)
        {
            double[] Y, X;

            if (pressure > bubblePointPressure)
            {
                Y = oilUnderSaturatedData[0]; X = oilUnderSaturatedData[1];
                if (pressure < Y[Y.Length - 1])
                {
                    return LookUp(Y, X, pressure);
                }
                else
                {
                    return Extrapolate(Y, X, pressure);
                }
            }
            else
            {
                Y = oilData[0]; X = oilData[1];
                return LookUp(Y, X, pressure);
            }

        }
        private double GetWaterFVF(double pressure)
        {
            double deltaP = pressure - waterData[0];
            double Bw = (2 * waterData[1] - waterData[2] * deltaP * waterData[1]) / (waterData[2] * deltaP + 2);
            return Bw;
        }
        private double GetGasFVF(double pressure)
        {
            double[] Y, X;
            Y = gasData[0]; X = gasData[1];

            double FVF = LookUp(Y, X, pressure);

            return FVF > 0 ? FVF : GetGasFVF(gasData[0].Last());
        }

        private double GetOilViscosity(double pressure)
        {
            double[] Y, X;

            if (pressure > bubblePointPressure)
            {
                Y = oilUnderSaturatedData[0]; X = oilUnderSaturatedData[2];
                return Extrapolate(Y, X, pressure);
            }
            else
            {
                Y = oilData[0]; X = oilData[2];
                return LookUp(Y, X, pressure);
            }
        }
        private double GetWaterViscosity(double pressure)
        {
            if (waterData[4] == 0)
            {
                return waterData[3];
            }
            else
            {
                throw new NotImplementedException();
            }
        }
        private double GetGasViscosity(double pressure)
        {
            double[] Y, X;
            Y = gasData[0]; X = gasData[2];

            double FVF = LookUp(Y, X, pressure);

            return pressure <= gasData[0].Last() ? FVF : Extrapolate(Y, X, pressure);
        }

        public double GetOilDensity(double pressure, double oil_fvf, double Rs)
        {
            double density = (surfaceDensities[0] + surfaceDensities[2] * Rs) / oil_fvf;
            return density;
        }
        public double GetWaterDensity(double pressure, double fvf)
        {
            double density = surfaceDensities[1] / fvf;
            return density;
        }
        public double GetGasDensity(double pressure, double fvf, double Rv)
        {
            double density = (surfaceDensities[2] + Rv * surfaceDensities[0]) / fvf;
            return density;
        }

        private double GetRso(double pressure)
        {
            double[] Y, X;

            if (pressure > bubblePointPressure)
            {
                var temp = oilUnderSaturatedData[3][0];
                return temp;
            }
            else
            {
                Y = oilData[0]; X = oilData[3];
                var temp = LookUp(Y, X, pressure);
                return temp;
            }

        }
        private double GetRsw(double pressure)
        {
            return 0;
        }


        // uses a Newton-Rapshon iteration to calculate the natural gas Z factor based on Dranchuk and Abu-Kassem equations
        // this method has no use currently. I developed it for an earlier test scenario and i see it's wasteful to delete it.
        private static double CalculateZ_Factor(double Pc, double Tc, double P, double T, double z_initial)
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

        public double GetAverageOilGravity(BaseBlock block, BaseBlock neighbor)
        {
            return 0.5 * (GetOilDensity(block.P[0], block.Bo[0], block.Rso[0]) + GetOilDensity(neighbor.P[0], neighbor.Bo[0], neighbor.Rso[0])) * Global.gamma_c * Global.alpha;
        }

        public double GetAverageGasGravity(BaseBlock block, BaseBlock neighbor)
        {
            return 0.5 * (GetGasDensity(block.P[0], block.Bg[0], block.Rvo[0]) + GetGasDensity(neighbor.P[0], neighbor.Bg[0], neighbor.Rvo[0])) * Global.gamma_c * Global.alpha;
        }

        public double GetAverageWaterGravity(BaseBlock block, BaseBlock neighbor)
        {
            return 0.5 * (GetWaterDensity(block.P[0], block.Bw[0]) + GetWaterDensity(neighbor.P[0], neighbor.Bw[0])) * Global.gamma_c * Global.alpha;
        }
    }
}
