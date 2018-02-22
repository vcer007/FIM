using FIM.Core;
using MathNet.Numerics.LinearAlgebra.Storage;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Mathematics
{
    struct SparseData
    {
        double value;
        int columnIndex;
        int rowIndex;
    }

    static class TempSparseMatrix
    {
        public static List<double> values;
        public static List<int> columnIndices;
        public static List<int> rowIndices;
        public static int counter;

        static TempSparseMatrix()
        {
            Reset();
        }

        public static void Reset()
        {
            counter = 0;
            values = new List<double>();
            columnIndices = new List<int>();
            rowIndices = new List<int>();
        }

        public static void SetAt(int row, int column, double value)
        {
            if (value != 0)
            {
                counter++;
                values.Add(value);
                rowIndices.Add(row);
                columnIndices.Add(column);
            }
        }

        public static List<double> GetValues()
        {
            return values;
        }

        public static List<int> GetColumnIndices()
        {
            return columnIndices;
        }

        public static List<int> GetRowIndices()
        {
            List<int> tempList = new List<int>();

            tempList.Add(0);
            for (int i = 1; i < values.Count; i++)
            {
                if (rowIndices[i] > rowIndices[i-1])
                {
                    tempList.Add(i);
                }
            }

            // add the final count
            tempList.Add(values.Count /*+ 1*/);

            return tempList;
        }

        public static void Sort()
        {
            bool swapped = true;
            for (int i = 0; i < values.Count; i++)
            {
                swapped = false;
                for (int j = 0; j < values.Count - 1 - i; j++)
                {
                    if (rowIndices[j] > rowIndices[j + 1] || (rowIndices[j] == rowIndices[j + 1] && columnIndices[j] > columnIndices[j + 1]))
                    {
                        Swap(j, rowIndices);
                        Swap(j, columnIndices);
                        Swap(j, values);
                        swapped = true;
                    }
                }
                if (!swapped)
                {
                    break;
                }
            }
        }

        static void Swap<T>(int firstIndex, List<T> array)
        {
            T temp = array[firstIndex];
            array[firstIndex] = array[firstIndex + 1];
            array[firstIndex + 1] = temp;
        }
    }

    static class MathNetSparseMatrixInitializer
    {
        static double[] _values;
        static int[] _columnIndices;
        static int[] _rowIndices;
        static int _numberOfNoneZeroElements;

        public static void Initialize(SimulationData data)
        {

            BaseBlock block;
            int total_count = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                total_count += 3 * 3;

                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] >= 0)
                    {
                        total_count += 3 * 3;
                    }
                }
            }

            _values = new double[total_count];
            _columnIndices = new int[total_count];
            _rowIndices = new int[total_count];
        }

        public static void SetAtIndex(int index, int row, int column, double value)
        {
            _values[index] = value;
            _rowIndices[index] = row;
            _columnIndices[index] = column;
        }

        public static List<double> GetValues()
        {
            var values_list = new List<double>();
            for (int i = 0; i < _values.Length; i++)
            {
                if (_values[i] != 0)
                {
                    values_list.Add(_values[i]);
                }
            }
            _numberOfNoneZeroElements = values_list.Count;
            return values_list;
        }

        public static List<int> GetRowPointers()
        {
            var row_pointers_list = new List<int>();
            row_pointers_list.Add(_rowIndices[0]);
            for (int i = 1; i < _values.Length; i++)
            {
                if (_values[i] != 0 && _rowIndices[i] != _rowIndices[i-1])
                {
                    row_pointers_list.Add(_rowIndices[i]);
                }
            }
            row_pointers_list.Add(_numberOfNoneZeroElements);

            return row_pointers_list;
        }

        public static List<int> GetColumnPointers()
        {
            var column_pointers_list = new List<int>();
            for (int i = 1; i < _values.Length; i++)
            {
                if (_values[i] != 0)
                {
                    column_pointers_list.Add(_columnIndices[i]);
                }
            }

            return column_pointers_list;
        }
    }

    static class MathNetSparseMatrixHelper
    {
        public static void InitializeStructure(this MathNet.Numerics.LinearAlgebra.Double.SparseMatrix Jacobi, SimulationData data, double value)
        {
            BaseBlock block;
            int counter = 0;
            int total_count = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                #region Oil
                // with respect to P
                Jacobi.At(counter, data.phases.Length * block.index, value);
                // with respect to Sg
                Jacobi.At(counter, data.phases.Length * block.index + 1, value);
                // with respect to Sw
                Jacobi.At(counter, data.phases.Length * block.index + 2, value);

                total_count += 3;

                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] >= 0)
                    {
                        // with respect to P
                        Jacobi.At(counter, data.phases.Length * block.neighborBlocksIndices[j], value);
                        // with respect to Sg
                        Jacobi.At(counter, data.phases.Length * block.neighborBlocksIndices[j] + 1, value);
                        // with respect to Sw
                        Jacobi.At(counter, data.phases.Length * block.neighborBlocksIndices[j] + 2, value);

                        total_count += 3;
                    }
                }
                #endregion
                #region Gas
                // with respect to P
                Jacobi.At(counter + 1, data.phases.Length * block.index, value);
                // with respect to Sg
                Jacobi.At(counter + 1, data.phases.Length * block.index + 1, value);
                // with respect to Sw
                Jacobi.At(counter + 1, data.phases.Length * block.index + 2, value);

                total_count += 3;

                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] >= 0)
                    {
                        // with respect to P
                        Jacobi.At(counter + 1, data.phases.Length * block.neighborBlocksIndices[j], value);
                        // with respect to Sg
                        Jacobi.At(counter + 1, data.phases.Length * block.neighborBlocksIndices[j] + 1, value);
                        // with respect to Sw
                        Jacobi.At(counter + 1, data.phases.Length * block.neighborBlocksIndices[j] + 2, value);

                        total_count += 3;
                    }
                }
                #endregion
                #region Water
                // with respect to P
                Jacobi.At(counter + 2, data.phases.Length * block.index, value);
                // with respect to Sg
                Jacobi.At(counter + 2, data.phases.Length * block.index + 1, value);
                // with respect to Sw
                Jacobi.At(counter + 2, data.phases.Length * block.index + 2, value);

                total_count += 3;

                for (int j = 0; j < block.neighborBlocksIndices.Length; j++)
                {
                    if (block.neighborBlocksIndices[j] >= 0)
                    {
                        // with respect to P
                        Jacobi.At(counter + 2, data.phases.Length * block.neighborBlocksIndices[j], value);
                        // with respect to Sg
                        Jacobi.At(counter + 2, data.phases.Length * block.neighborBlocksIndices[j] + 1, value);
                        // with respect to Sw
                        Jacobi.At(counter + 2, data.phases.Length * block.neighborBlocksIndices[j] + 2, value);

                        total_count += 3;
                    }
                }
                #endregion

                counter += data.phases.Length;
            }

            ((SparseCompressedRowMatrixStorage<double>)Jacobi.Storage).RowPointers[data.grid.Length * data.phases.Length] = total_count;

        }

        public static void SetValueUchecked(this MathNet.Numerics.LinearAlgebra.Double.SparseMatrix Jacobi, int position, double value)
        {
            ((SparseCompressedRowMatrixStorage<double>)Jacobi.Storage).Values[position] = value;
        }
    }
}
