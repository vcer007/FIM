using FIM.Core;
using FIM.MaterialBalance;
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
        /// The indices, in the <see cref="SimulationData.grid"/>, of the blocks whose properties are to be reported.
        /// </summary>
        public List<string[]> blocksXYZIndices = new List<string[]>();

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
        public void Write(double currentTime, SimulationData data, bool headLine = false, int counter = 0, int milli_seconds = 0)
        {
            StringBuilder output = new StringBuilder();

            if (headLine)
            {
                output.Append("Time".PadRight(Global.padding));

                for (int i = 0; i < wellKeyWords.Count; i++)
                {
                    for (int j = 0; j < wellsIndices[i].Length; j++)
                    {
                        string temp = (wellKeyWords[i]).PadRight(Global.padding);
                        output.Append(temp);
                    }
                }

                for (int i = 0; i < blockKeyWords.Count; i++)
                {
                    //for (int j = 0; j < blocksIndices[i].Length; j++) output.Append( (blockKeyWords[i] + "_BLOCK@" + blocksIndices[i][j]).PadRight(Global.padding));
                    for (int j = 0; j < blocksIndices[i].Length; j++) output.Append((blockKeyWords[i]).PadRight(Global.padding));
                }

                foreach (string keyword in singleKeyWords)
                {
                    output.Append(keyword.PadRight(Global.padding));
                }

                output.AppendLine();
                output.AppendLine();

                output.Append(string.Empty.PadRight(Global.padding));

                for (int i = 0; i < wellKeyWords.Count; i++)
                {
                    for (int j = 0; j < wellsIndices[i].Length; j++)
                    {
                        string temp = data.wells.First(x => x.index == wellsIndices[i][j]).name.PadRight(Global.padding);
                        output.Append(temp);
                    }
                }

                output.AppendLine();
                output.Append(string.Empty.PadRight(Global.padding));
                for (int i = 0; i < wellKeyWords.Count; i++)
                {
                    for (int j = 0; j < wellsIndices[i].Length; j++)
                    {
                        output.Append(string.Empty.PadRight(Global.padding));
                    }
                }

                for (int i = 0; i < blockKeyWords.Count; i++)
                {
                    for (int j = 0; j < blocksIndices[i].Length; j++) output.Append((string.Join("  ", blocksXYZIndices[i][j].Split(','))).PadRight(Global.padding));
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
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells.First(w => w.index == wellsIndices[i][j]).BHP[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;
                        
                        case "WOPR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells.First(w => w.index == wellsIndices[i][j]).GetTotalOilProduction().ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WOPRF":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells.First(w => w.index == wellsIndices[i][j]).q_oil[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WOPRS":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells.First(w => w.index == wellsIndices[i][j]).q_vap_oil[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WWPR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append(data.wells.First(w => w.index == wellsIndices[i][j]).q_water[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGPRF":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells.First(w => w.index == wellsIndices[i][j]).q_free_gas[0] * Global.a / 1000).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGPRS":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells.First(w => w.index == wellsIndices[i][j]).q_solution_gas[0] * Global.a / 1000).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGOR":
                            //for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells.First(w => w.index == wellsIndices[i][j]).q_solution_gas[0] / data.wells.First(w => w.index == wellsIndices[i][j]).q_oil[0]).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells.First(w => w.index == wellsIndices[i][j]).GetProducingGOR()).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "WGPR":
                            for (int j = 0; j < wellsIndices[i].Length; j++) output.Append((data.wells.First(w => w.index == wellsIndices[i][j]).q_solution_gas[0] * Global.a / 1000 + data.wells.First(w => w.index == wellsIndices[i][j]).q_free_gas[0] * Global.a / 1000).ToString(Global.decimalPlaces).PadRight(Global.padding));
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

                        case "BOVIS":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append(data.grid[blocksIndices[i][j]].viscosityOil[0].ToString(Global.decimalPlaces).PadRight(Global.padding));
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

                        case "BODEN":
                            for (int j = 0; j < blocksIndices[i].Length; j++) { BaseBlock block = data.grid[blocksIndices[i][j]]; output.Append((data.pvt.GetOilDensity(block.P[0], block.Bo[0], block.Rso[0])).ToString(Global.decimalPlaces).PadRight(Global.padding)); };
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
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append((data.grid[blocksIndices[i][j]].Rso[0]*Global.a/1000).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "BRV":
                            for (int j = 0; j < blocksIndices[i].Length; j++) output.Append((data.grid[blocksIndices[i][j]].Rvo[0] * 1000 / Global.a).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        default:
                            break;
                    }
                }

                foreach (string keyword in singleKeyWords)
                {

                    switch (keyword)
                    {
                        case "FOPR":
                            output.Append(data.wells.Sum(w => w.q_oil[0]).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "MBEO":
                            output.Append((data.MBE_Oil / data.MBE_Tolerance).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "MBEW":
                            output.Append((data.MBE_Water / data.MBE_Tolerance).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "MBEG":
                            output.Append((data.MBE_Gas / data.MBE_Tolerance).ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "NEWTON":
                            output.Append(counter.ToString().PadRight(Global.padding));
                            break;

                        case "TCPUTS":
                            output.Append((milli_seconds / 1000.0).ToString().PadRight(Global.padding));
                            break;

                        case "FGIP":
                            var temp = MBE.GetGIP(data, 0);
                            output.Append(temp.ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "FGIPG":
                            temp = MBE.FreeGasIP(data, 0);
                            output.Append(temp.ToString(Global.decimalPlaces).PadRight(Global.padding));
                            break;

                        case "FGIPL":
                            temp = MBE.SolubleGasIP(data, 0);
                            output.Append(temp.ToString(Global.decimalPlaces).PadRight(Global.padding));
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
                if (!console)
                {
                    if (headLine)
                    {
                        Console.WriteLine("Added the headlines to the output file.");
                    }
                    else
                    {
                        string numberOrder = string.Empty;
                        string stepsCounterAsString = stepsCounter.ToString();
                        string lastTwo = stepsCounter > 10 ? stepsCounterAsString.Substring(stepsCounterAsString.Length - 2, 2) : string.Empty;
                        char lastOne = stepsCounterAsString.Last();
                        if ((lastTwo == "11") || (lastTwo == "12") || (lastTwo == "13"))
                            numberOrder = "th";
                        else
                        {
                            numberOrder = lastOne == '3' ? "rd" : lastOne == '2' ? "nd" : lastOne == '1' ? "st" : "th";
                        }

                        Console.WriteLine("Added the {0}{1} time step output to the file.", stepsCounter, numberOrder);
                        stepsCounter++;
                    }
                }
                

                using (StreamWriter fileWriter = new StreamWriter(data.title + ".OUT", !headLine))
                {
                    fileWriter.WriteLine(output);
                }
            }
        }
    }
}
