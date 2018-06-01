using System;
using FIM.Core;

namespace FIM.FluidData
{
    /// <summary>
    /// The Special Core Analysis data storage class.
    /// </summary>
    /// <remarks>
    /// <para>This class acts as a container to store all input PVT data.</para>
    /// <para>It also contains methods for interpolation to get PVT data at different values of pressure.</para>
    /// </remarks>
    /// <seealso cref="PVT"/>
    /// <seealso cref="RockData.PorosityCalculator"/>
    public class SCAL
    {
        // this array is used to store the input relative permeability data.
        private double[][] kr_Data;
        private double[][] sgfn;
        private double[][] swfn;
        private double[][] sof3;

        /// <summary>
        /// Initializes a new instance of the <see cref="SCAL"/> class.
        /// </summary>
        /// <remarks>
        /// This is used only once to get the Kr data stored within the Kr class instance.
        /// </remarks>
        /// <param name="kr_Data">
        /// takes a matrix (double[][] array) as an input 
        /// where the first column is S then Krg, Kro then Krw.
        /// </param>
        public SCAL(double[][] kr_Data)
        {
            this.kr_Data = kr_Data;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="SCAL"/> class.
        /// </summary>
        /// <remarks>
        /// This is an empty constructor of a new SCAL instance.
        /// The <see cref="Initialize(double[][])"/> method must be used to assign the data.
        /// </remarks>
        public SCAL()
        {

        }

        /// <summary>
        /// Initializes the specified Kr data.
        /// </summary>
        /// <remarks>
        /// This method is only used when a <see cref="SCAL"/> instance is created using the empty constructor.
        /// The empty constructor is used to create a <see cref="SCAL"/> object "instance" without specifying its data.
        /// This should be taken with so much care.
        /// </remarks>
        /// <param name="kr_Data">The Kr data.</param>
        public void Initialize(double[][] kr_Data)
        {
            this.kr_Data = kr_Data;
        }

        /// <summary>
        /// Gets the relative permeability for a certain phase.
        /// </summary>
        /// <param name="phase">The phase.</param>
        /// <param name="saturation">The saturation.</param>
        /// <returns>The value of the relative permeability</returns>
        /// <seealso cref="Global.Phase"/>
        //public double GetKr(Global.Phase phase, double saturation)
        //{
        //    switch (phase)
        //    {
        //        case Global.Phase.Water:
        //            return GetKrw(saturation);
        //        case Global.Phase.Oil:
        //            return GetKro(saturation);
        //        case Global.Phase.Gas:
        //            return GetKrg(saturation);
        //        default:
        //            return 1;
        //    }
        //}


        //Internal helper method to do the table lookups and interpolations
        //Method Name: lookUp
        //Objectives: this is a private (internal method for the internal use of this class only) method that makes table lookup interpolation
        //Inputs: two variables representing the X and Y arrays "The two columns of the table" and the Y value
        //Outputs: a variable that represents the X value corresponding to the particular Y value
        // this function assumes that array elements are sorted ascendingly
        private double LookUp(double[] data_y, double[] data_x, double y)
        {
            double temp = data_x[data_x.Length - 1];

            if (y > data_y[data_y.Length - 1])
            {
                return temp;
            }

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
                    return x;
                }
            }

            return temp;
        }


        //Method Name: getKr "where Kr may be Kro, Krw or Krg"
        //Objectives: calculate the corresponding Kr value at a certain saturation
        //Inputs: a variable representing the value of the saturation
        //Outputs: the corresponding Kr value

        public double GetKro(double sg, double sw, double swco)
        {
            double so = 1 - sg - sw;
            so = so < 0 ? 0 : so;
            double krog = LookUp(sof3[0], sof3[2], so);
            double krow = LookUp(sof3[0], sof3[1], so);

            double temp = (sg + sw - swco);
            double kro = temp == 0 ? 1 : (sg * krog + (sw - swco) * krow) / temp;
            return kro > 1 ? 1 : kro;
        }

        public double GetKrg(double saturation)
        {
            double kr = LookUp(sgfn[0], sgfn[1], saturation);
            return kr > 1 ? 1 : kr;
        }

        public double GetKrw(double saturation)
        {
            //return 0;
            if (saturation < swfn[0][0])
            {
                return 0;
            }
            double kr = LookUp(swfn[0], swfn[1], saturation);
            return kr > 1 ? 1 : kr;
        }

        public double GetWaterCapillaryPressure(double saturation)
        {
            if (saturation < swfn[0][0])
            {
                return swfn[2][0];
            }

            return LookUp(this.swfn[0]/*water saturation*/, this.swfn[2]/*capillary*/, saturation);
        }

        public double GetGasCapillaryPressure(double saturation)
        {
            return LookUp(this.sgfn[0]/*water saturation*/, this.sgfn[2]/*capillary*/, saturation);
        }

        internal void Initialize(double[][] sgfn, double[][] swfn, double[][] sof3)
        {
            this.sgfn = sgfn;
            this.swfn = swfn;
            this.sof3 = sof3;
        }
    }
}
