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


        //public double GetKro(double sg, double sw, double swco)
        //{
        //    double[] Y, X;
        //    Y = kr_Data[0]; X = kr_Data[2];

        //    return LookUp(Y, X, sg);
        //}

        public double GetKro(double sg, double sw, double swco)
        {
            double[] Y, X;
            Y = kr_Data[0]; X = kr_Data[2];

            Y = new double[] { 0, 0.18, 0.28, 0.38, 0.43, 0.48, 0.58, 0.63, 0.68, 0.76, 0.83, 0.86, 0.879, 0.88 };
            X = new double[] { 0, 0, 0.0001, 0.001, 0.01, 0.021, 0.09, 0.2, 0.35, 0.7, 0.98, 0.997, 1, 1 };

            double so = 1 - sg - sw;
            so = so < 0 ? 0 : so;
            double krog = LookUp(Y, X, so);
            double krow = LookUp(Y, X, so);

            //return LookUp(Y, X, saturation);
            double temp = (sg + sw - swco);
            double kro = temp == 0 ? 1 : (sg * krog + (sw - swco) * krow) / temp;
            return kro > 1 ? 1 : kro;
        }

        public double GetKrg(double saturation)
        {
            double[] Y, X;
            Y = kr_Data[0]; X = kr_Data[1];

            double kr = LookUp(Y, X, saturation);
            return kr > 1 ? 1 : kr;
        }

        public double GetKrw(double saturation)
        {

            double[] Y, X;
            //Y = kr_Data[0]; X = kr_Data[3];

            Y = new double[] {0, 0.12, 1};
            X = new double[] {0, 0, 0.00001/*1*/ };

            double kr = LookUp(Y, X, saturation);
            return 0;
            //return kr > 1 ? 1 : kr;
        }
    }
}
