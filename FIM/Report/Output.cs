using FIM.Core;
using FIM.Well;

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Report
{
    /// <summary>
    /// This class contains methods for handling simulation run output to files or console.
    /// </summary>
    public class Output
    {
        // stores the number of steps simulated.
        int stepsCounter = 1;

        /// <summary>
        /// Whether or not the <see cref="Write"/> method should store reports to an output file.
        /// </summary>
        public bool file = false;

        /// <summary>
        /// Whether or not the <see cref="Write"/> method should display reports in the <see cref="Console"/>
        /// </summary>
        public bool console = false;

        /// <summary>
        /// The indices, in the <see cref="SimulationData.grid"/>, of the blocks whose properties are to be reported.
        /// </summary>
        public List<int[]> blocksIndices = new List<int[]>();

        /// <summary>
        /// A list of the grid block properties to be reported.
        /// </summary>
        public List<string> blockKeyWords= new List<string>();

        /// <summary>
        /// The indices, in the <see cref="SimulationData.wells"/>, of the wells whose properties are to be reported.
        /// </summary>
        public List<int[]> wellsIndices = new List<int[]>();

        /// <summary>
        /// A list of the grid block properties to be reported.
        /// </summary>
        public List<string> wellKeyWords = new List<string>();

        /// <summary>
        /// A list of single key words to be reported.
        /// </summary>
        /// <remarks>
        /// This list includes keywords that have no further data required "indices" like the material balance errors.
        /// </remarks>
        public List<string> singleKeyWords = new List<string>();

        /// <summary>
        /// Writes the specified reporting data to the specified locations.
        /// </summary>
        /// <param name="currentTime">The current time of the simulation run.</param>
        /// <param name="data">The <see cref="SimulationData"/> instance containing access to all the run data.</param>
        /// <param name="headLine">Indicates if this is the first line to be reported. This is used to write the headlines of the report columns.</param>
        public void Write(double currentTime, SimulationData data, bool headLine = false)
        {
            StringBuilder output = new StringBuilder();

            if (headLine)
            {
                output.Append("Time".PadRight(Global.padding));

                for (int i = 0; i < wellKeyWords.Count; i++)
                {
                    for (int j = 0; j < wellsIndices[i].Length; j++)
                    {
                        string temp = (wellKeyWords[i] + "_" + data.wells.Where(x => x.index == wellsIndices[i][j]).ToArray()[0].name).PadRight(Global.padding);
                        output.Append(temp);
                    }
                }

                for (int i = 0; i < blockKeyWords.Count; i++)
                {
                    for (int j = 0; j < blocksIndices[i].Length; j++) output.Append( (blockKeyWords[i] + "_BLOCK@" + blocksIndices[i][j]).PadRight(Global.padding));
                }

                foreach (string keyword in singleKeyWords)
                {
                    output.Append(keyword.PadRight(Global.padding));
                }
            }
            else
            {
                output.Append(currentTime.ToString().PadRight(Global.padding));

                // well
                for (int i = 0; i < wellKeyWords.Count; i++)
                {
                    switch (wellKeyWords[i])
                    {
                        case "WBHP":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells[wellsIndices[i][j]].BHP[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WOPR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells[wellsIndices[i][j]].q_oil[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WWPR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells[wellsIndices[i][j]].q_water[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGPRF":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells[wellsIndices[i][j]].q_free_gas[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGPRS":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells[wellsIndices[i][j]].q_solution_gas[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGOR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells[wellsIndices[i][j]].q_solution_gas[0] / data.wells[j].q_oil[0]).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGPR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells[wellsIndices[i][j]].q_solution_gas[0] + data.wells[j].q_free_gas[0]).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        default:
                            break;
                    }
                }
                
                // block
                for (int i = 0; i < blockKeyWords.Count; i++)
                {
                    switch (blockKeyWords[i])
                    {
                        case "BOFVF":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Bo[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BWFVF":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Bw[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BGFVF":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Bg[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BKRO":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Kro[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BKRW":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Krw[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BKRG":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Krg[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BPR":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].P[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BOSAT":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].So[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BWSAT":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Sw[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BGSAT":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Sg[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BRS":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].Rso[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        default:
                            break;
                    }
                }

                foreach (string keyword in singleKeyWords)
                {
                    switch (keyword)
                    {
                        case "MBEO":
                            output.Append(data.MBE_Oil.ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "MBEW":
                            output.Append(data.MBE_Water.ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "MBEG":
                            output.Append(data.MBE_Gas.ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        default:
                            break;
                    }
                }
            }

            if (console)
            {
                Console.WriteLine(output);
            }

            if (file)
            {
                if (headLine)
                {
                    Console.WriteLine("Added the headlines to the output file.");
                }
                else
                {
                    string numberOrder = stepsCounter == 1 ? "st" : stepsCounter == 2 ? "nd" : stepsCounter == 3 ? "rd" : "th";
                    Console.WriteLine("Added the {0}{1} time step output to the file.", stepsCounter, numberOrder);
                    stepsCounter++;
                }
                

                using (StreamWriter fileWriter = new StreamWriter(data.title + ".OUT", !headLine))
                {
                    fileWriter.WriteLine(output);
                }
            }
        }
    }
}
