using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace FIM.Parser
{
    /// <summary>
    /// This class contains helper methods useful for parsing input text files.
    /// </summary>
    public static class Helper
    {
        /// <summary>
        /// Determines whether the line contains any of the specified key words.
        /// </summary>
        /// <param name="line">The line.</param>
        /// <param name="KeyWords">The key words.</param>
        /// <returns></returns>
        public static bool ContainsAnyOf(this string line, string[] KeyWords)
        {
            for (int i = 0; i < KeyWords.Length; i++)
            {
                //if (line.Contains(KeyWords[i]))
                //{
                //    return true;
                //}

                if (Regex.Match(line, "\\b" + KeyWords[i] + "\\b").Success)
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Reads the file specified at the path.
        /// </summary>
        /// <param name="filePath">The file path.</param>
        /// <returns>
        /// A list of all the lines in the file.
        /// </returns>
        public static List<string> ReadFile(string filePath)
        {
            List<string> lines = new List<string>();

            // read all the lines.
            string[] file = File.ReadAllLines(filePath);

            // remove commented-out sections from the lines.
            for (int i = 0; i < file.Length; i++)
            {
                file[i] = Regex.Replace(file[i], @"\-(-+)\s?.*", "");
            }

            // remove back-slashes
            for (int i = 0; i < file.Length; i++)
            {
                file[i] = Regex.Replace(file[i], @"/", String.Empty);
            }

            // add non-empty lines to the list.
            foreach (string line in file)
            {
                // not empty.
                if (Regex.Match(line, @"[A-Za-z0-9]+").Success)
                {
                    // add to the lines list
                    lines.Add(line);
                }
            }

            // remove extra spaces at the beginning of the line
            // remove extra characters after the main data file sections

            return lines;
        }


        /// <summary>
        /// Gets the data at a specified line.
        /// </summary>
        /// <param name="line">The line.</param>
        /// <param name="size">The size of the data in the line.</param>
        /// <returns></returns>
        public static double[] GetData(string line, int size = 1)
        {
            string[] data = Regex.Split(line, @"\b\s+\b");
            double[] numericData;

            if (data.Length == 1)
            {
                numericData = Enumerable.Repeat(double.Parse(data[0]), size).ToArray();
            }
            else
            {
                numericData = new double[data.Length];
                for (int i = 0; i < numericData.Length; i++)
                {
                    numericData[i] = double.Parse(data[i]);
                }
            }

            return numericData;
        }

        /// <summary>
        /// Looks for a dfinitive keyword in the whole section and returns its associated data.
        /// </summary>
        /// <param name="keyWord">The key word.</param>
        /// <param name="section">The section.</param>
        /// <param name="size">The size of the data in the line.</param>
        /// <returns></returns>
        public static double[] GetData(string keyWord, List<string> section, int size = 1)
        {
            int index = section.FindIndex(x => x.Contains(keyWord));

            string[] data = Regex.Split(section[index + 1], @"\b\s+\b");
            double[] numericData;

            if (data.Length == 1)
            {
                numericData = Enumerable.Repeat(double.Parse(data[0]), size).ToArray();
            }
            else
            {
                numericData = new double[data.Length];
                for (int i = 0; i < numericData.Length; i++)
                {
                    numericData[i] = double.Parse(data[i]);
                }
            }

            return numericData;
        }

        /// <summary>
        /// Gets the data as string.
        /// </summary>
        /// <remarks>
        /// This method does not convert the extracted data into a numerical value.
        /// This is useful when the extracted data itself contains keywords.
        /// </remarks>
        /// <param name="line">The line.</param>
        /// <param name="size">The size of the data in the line.</param>
        /// <returns></returns>
        public static string[] GetDataAsString(string line, int size = 1)
        {
            string[] data = Regex.Split(line, @"\s+").Where(x => !String.IsNullOrEmpty(x)).ToArray();

            return data;
        }

        /// <summary>
        /// Gets the indices "x, y, z" in multiple lines.
        /// </summary>
        /// <param name="keyWord">The key word to look for.</param>
        /// <param name="allKeyWords">All key words to detect the end of the keyword.</param>
        /// <param name="section">All the lines in the section.</param>
        /// <returns>An array in which each element represents a line converted to a 3 elements array of the indices "x, y, z".</returns>
        public static int[][] GetMultipleRowsIndices(string keyWord, string[] allKeyWords,List<string> section)
        {
            int index = section.FindIndex(x => x.Contains(keyWord)) + 1;

            int index_end = index;
            for (int i = index; i < section.Count; i++)
            {
                if (section[i].ContainsAnyOf(allKeyWords))
                {
                    index_end = i - 1;
                    break;
                }
                else if (i == section.Count - 1)
                {
                    index_end = i;
                    break;
                }
            }

            int[][] data = new int[index_end - index + 1][];
            //for (int i = index; i <= index_end; i++)
            //{
            //    data[i] = GetData(section[index]).Select(x => (int)x).ToArray();
            //}
            for (int i = 0; i < data.Length; i++)
            {
                //data[i] = GetData(section[index]).Select(x => (int)x).ToArray();
                data[i] = Regex.Split(section[index], @"\s+").Where(x => !String.IsNullOrEmpty(x)).Select(x => int.Parse(x) - 1).ToArray();
                index++;
            }

            return data;
        }
    }
}
