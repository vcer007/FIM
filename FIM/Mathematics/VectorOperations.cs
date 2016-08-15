using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Mathematics
{
    class VectorOperations
    {
        public static double[] Subtract(double[] A, double[] B)
        {
            double[] temp = new double[A.Length];
            for (int i = 0; i < temp.Length; i++)
            {
                temp[i] = A[i] - B[i];
            }
            return temp;
        }

        public static double[] Add(double[] A, double[] B)
        {
            double[] temp = new double[A.Length];
            for (int i = 0; i < temp.Length; i++)
            {
                temp[i] = A[i] + B[i];
            }
            return temp;
        }

        public static double Dot(double[] A, double[] B)
        {
            double temp = 0;
            for (int i = 0; i < A.Length; i++)
            {
                temp += A[i] * B[i];
            }
            return temp;
        }

        public static double[] Multiply(double a, double[] B)
        {
            double[] temp = new double[B.Length];
            for (int i = 0; i < temp.Length; i++)
            {
                temp[i] = a * B[i];
            }
            return temp;
        }

        public static double[] SumMultiply(List<double> b, List<double[]> B)
        {
            double[] temp = new double[B[0].Length];
            for (int i = 0; i < b.Count; i++)
            {
                temp = Add(temp, Multiply(b[i], B[i]));
            }
            return temp;
        }

    }
}
