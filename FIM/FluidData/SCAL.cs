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
        /// Gets the relative permeability for a certain phase.
        /// </summary>
        /// <param name="phase">The phase.</param>
        /// <param name="saturation">The saturation.</param>
        /// <returns>The value of the relative permeability</returns>
        /// <seealso cref="Global.Phase"/>
        public double GetKr(Global.Phase phase, double saturation)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return GetKrw(saturation);
                case Global.Phase.Oil:
                    return GetKro(saturation);
                case Global.Phase.Gas:
                    return GetKrg(saturation);
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


        //Method Name: getKr "where Kr may be Kro, Krw or Krg"
        //Objectives: calculate the corresponding Kr value at a certain saturation
        //Inputs: a variable representing the value of the saturation
        //Outputs: the corresponding Kr value

        private double GetKro(double saturation)
        {
            double[] Y, X;
            Y = kr_Data[0]; X = kr_Data[2];

            return LookUp(Y, X, saturation);
        }

        private double GetKrg(double saturation)
        {
            double[] Y, X;
            Y = kr_Data[0]; X = kr_Data[1];

            return LookUp(Y, X, saturation);
        }

        private double GetKrw(double saturation)
        {
            double[] Y, X;
            Y = kr_Data[0]; X = kr_Data[3];

            return LookUp(Y, X, saturation);
        }
    }
}
