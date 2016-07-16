using FIM.Core;

namespace FIM.Fluid
{
    //Method_1 Name: lookUp
    //Objectives: this is a private (internal method for the internal use of this class only) method that makes table lookup interpolation
    //Inputs: two variables representing the X and Y arrays "The two columns of the table" and the Y value
    //Outputs: a variable that represents the X value corresponding to the particular Y value

    //Method_2 Name: getKr "where Kr may be Kro, Krw or Krg"
    //Objectives: calculate the corresponding Kr value at a certain saturation
    //Inputs: a variable representing the value of the saturation
    //Outputs: the corresponding Kr value
    public class Kr
    {
        private double[][] Kr_data;

        //Objectives: Public constructor of the class. This is used only once to get the Kr data stored within the Kr class instance
        //Inputs: takes a matrix (double[][] array) as an input where the first column is S then Krg, Kro then Krw
        public Kr(double[][] kr_data)
        {
            this.Kr_data = kr_data;
        }


        public double getKr(Global.Phase phase, double saturation)
        {
            switch (phase)
            {
                case Global.Phase.Water:
                    return getKrw(saturation);
                case Global.Phase.Oil:
                    return getKro(saturation);
                case Global.Phase.Gas:
                    return getKrg(saturation);
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

        private double getKro(double saturation)
        {
            double[] Y, X;
            Y = Kr_data[0]; X = Kr_data[2];

            return lookUp(Y, X, saturation);
        }

        private double getKrg(double saturation)
        {
            double[] Y, X;
            Y = Kr_data[0]; X = Kr_data[1];

            return lookUp(Y, X, saturation);
        }

        private double getKrw(double saturation)
        {
            double[] Y, X;
            Y = Kr_data[0]; X = Kr_data[3];

            return lookUp(Y, X, saturation);
        }
    }
}
