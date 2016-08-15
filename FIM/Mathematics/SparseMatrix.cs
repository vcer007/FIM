using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Mathematics
{
    public class SparseMatrix
    {
        public double[] values;
        public int[] indicesX;
        public int[] indicesY;
        public int[] indicesDiagonal;

        public void SetElement(double value, int x, int y)
        {
            
        }

        public double[] Multiply(double[] vector)
        {
            double[] temp = new double[vector.Length];
            int counter = 0;

            for (int i = 0; i < temp.Length; i++)
            {
                // reset to zero.
                temp[i] = 0;
                while (counter < indicesY.Length && indicesY[counter] == i)
                {
                    temp[i] += values[counter] * vector[indicesX[counter]];
                    counter++;
                }
            }

            return temp;
        }

        public static SparseMatrix FromJagged(double[][] jaggedArray)
        {
            SparseMatrix matrix = new SparseMatrix();

            List<double> values = new List<double>();
            List<int> indicesX = new List<int>();
            List<int> indicesY = new List<int>();
            List<int> indicesDiagonal = new List<int>();

            for (int i = 0; i < jaggedArray.Length; i++)
            {
                for (int j = 0; j < jaggedArray[i].Length; j++)
                {
                    if (jaggedArray[i][j] != 0)
                    {
                        values.Add(jaggedArray[i][j]);
                        indicesX.Add(j);
                        indicesY.Add(i);
                    }

                    if (i == j)
                    {
                        indicesDiagonal.Add(values.Count - 1);
                    }
                }
            }

            matrix.values = values.ToArray();
            matrix.indicesX = indicesX.ToArray();
            matrix.indicesY = indicesY.ToArray();
            matrix.indicesDiagonal = indicesDiagonal.ToArray();

            return matrix;
        }
    }
}
