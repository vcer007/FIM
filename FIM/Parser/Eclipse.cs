using FIM.Core;
using FIM.FluidData;
using FIM.Initialization;
using FIM.RockData;
using FIM.Extensions;

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using FIM.Report;
using System.Text;

/// <summary>
/// This namespace contains classes for handling text files parsing.
/// </summary>
namespace FIM.Parser
{
    /// <summary>
    /// This class contains methods for reading input data files.
    /// </summary>
    /// <remarks>
    /// Keywords are kept in accordance with eclipse as possiblle as it could be.
    /// </remarks>
    public class Eclipse
    {
        // The data file sections keywords.
        static string[] dataFileSections = new string[] { "RUNSPEC", "GRID", "PROPS", "SOLUTION", "SUMMARY", "SCHEDULE" , "END"};

        // These lists are used to store data file sections eparately.
        // this allows for handling them separately.
        static List<string> RUNSPEC_section = new List<string>();
        static List<string> GRID_section = new List<string>();
        static List<string> PROPS_section = new List<string>();
        static List<string> SOLUTION_section = new List<string>();
        static List<string> SUMMARY_section = new List<string>();
        static List<string> SCHEDULE_section = new List<string>();

        // internal indexing variables for looking up different keywords in the file.
        static int startIndex = 0, endIndex = 0;
        static bool startIndexSet;

        // this variable is the return result of this class.
        // an initialized SimulationData instance.
        static SimulationData data = new SimulationData();


        static PorosityCalculator porosityCalculator = new PorosityCalculator();
        static PVT pvt = new PVT();
        static SCAL scal = new SCAL();
        static Output output = new Output();

        static double So, Sw, Sg, P, Swc;

        /// <summary>
        /// Reads the input files from the specified paths.
        /// </summary>
        /// <param name="dataFilePath">The data file path.</param>
        /// <param name="initializationFilePath">The initialization file path.</param>
        /// <remarks>
        /// If the files are located within the same directory "folder" as the simulator executable ".exe", using the names only is sufficient.
        /// </remarks>
        /// <returns>
        /// An initialized instance of the <see cref="SimulationData"/> class containing the data read from the files.
        /// </returns>
        public static SimulationData ReadInputData(string dataFilePath, string initializationFilePath)
        {
            // internal variable to store the file's lines as a list.
            List<string> lines;

            #region Initialization File

            lines = Helper.ReadFile(initializationFilePath);

            InitializationFile(lines);

            #endregion

            #region Data File

            lines = Helper.ReadFile(dataFilePath);

            // add the sections
            AddSection("RUNSPEC", RUNSPEC_section, lines, dataFileSections);
            AddSection("GRID", GRID_section, lines, dataFileSections);
            AddSection("PROPS", PROPS_section, lines, dataFileSections);
            AddSection("SOLUTION", SOLUTION_section, lines, dataFileSections);
            AddSection("SUMMARY", SUMMARY_section, lines, dataFileSections);
            AddSection("SCHEDULE", SCHEDULE_section, lines, dataFileSections);

            InitializeRUNSPEC(data, RUNSPEC_section);
            InitializeGRID(data, GRID_section, porosityCalculator);
            InitializePROPS(PROPS_section, pvt, scal, porosityCalculator);
            InitializeSOLUTION(SOLUTION_section);
            InitializeSCHEDULE(data, SCHEDULE_section);
            InitializeSUMMARY(SUMMARY_section);


            data.porosityCalculator = porosityCalculator;
            data.pvt = pvt;
            data.scal = scal;
            data.output = output;

            // update block properties
            for (int i = 0; i < data.grid.Length; i++)
            {
                data.grid[i].UpdateProperties(data, P, Sw, Sg, So, 0);
            }

            #endregion

            return data;            
        }


        // get data from input data file.

        private static void InitializeRUNSPEC(SimulationData data, List<string> section)
        {
            int index;
            // first : single-line keywords

            // phases
            List<Global.Phase> phases = new List<Global.Phase>();
            if (section.Contains("OIL")) phases.Add(Global.Phase.Oil);
            if (section.Contains("WATER")) phases.Add(Global.Phase.Water);
            if (section.Contains("GAS")) phases.Add(Global.Phase.Gas);
            if (section.Contains("DISGAS"))
            {
                data.solubleGasPresent = true;
            }
            if (section.Contains("VAPOIL"))
            {
                data.vaporizedOilPresent = true;
            }

            if (section.Contains("GRAVITY")) data.Gravity = true;

            data.phases = phases.ToArray();

            if (section.Contains("IMPES"))
            {
                data.solutionProcedure = Global.SolutionProcedure.IMPES;
            }
            else
            {
                data.solutionProcedure = Global.SolutionProcedure.FullyImplicit;
            }
            

            // second : two-lines keywords

            // title
            index = section.FindIndex(x => x.Contains("TITLE"));
            data.title = section[index + 1];

            // grid dimensions
            index = section.FindIndex(x => x.Contains("DIMENS"));
            string[] dimens = section[index + 1].Split(new char[] { }, StringSplitOptions.RemoveEmptyEntries);
            data.x = int.Parse(dimens[0]);
            data.y = int.Parse(dimens[1]);
            data.z = int.Parse(dimens[2]);

            // well dimensions
            index = section.FindIndex(x => x.Contains("WELLDIMS"));
            dimens = section[index + 1].Split(new char[] { }, StringSplitOptions.RemoveEmptyEntries);
            data.wells = new Well.BaseWell[int.Parse(dimens[0])];
        }

        private static void InitializeGRID(SimulationData data, List<string> section, PorosityCalculator porosityCalculator)
        {
            double[] x_dimensionLengths = new double[data.x];
            double[] y_dimensionLengths = new double[data.y];
            double[] heights = new double[data.z];
            double[] transmissibilitiyMultipliers = new double[data.z];
            double[] layersPermeabilities = new double[data.z];

            // grid blocks dimensions

            // x
            int index = section.FindIndex(x => x.Contains("DX"));
            x_dimensionLengths = Helper.GetData(section[index + 1], data.x);
            // y
            index = section.FindIndex(x => x.Contains("DY"));
            y_dimensionLengths = Helper.GetData(section[index + 1], data.y);
            // z
            index = section.FindIndex(x => x.Contains("DZ"));
            heights = Helper.GetData(section[index + 1], data.z);

            // permeabilities

            // permeability in the x direction
            index = section.FindIndex(x => x.Contains("PERMX"));
            layersPermeabilities = Helper.GetData(section[index + 1], data.z);

            // transmissibilities multipliers in the z direction only.
            index = section.FindIndex(x => x.Contains("MULTZ"));
            if (index != -1)
            {
                transmissibilitiyMultipliers = Helper.GetData(section[index + 1], data.z);
            }
            else
            {
                transmissibilitiyMultipliers = Enumerable.Repeat(1.0, data.z).ToArray();
            }

            // porosity
            index = section.FindIndex(x => x.Contains("PORO"));
            // only one value of porosity is initially assigned to the whole model.
            porosityCalculator.porosity_ref = Helper.GetData(section[index + 1], 1)[0];



            CartesianGrid cartesianGrid = new CartesianGrid(data.x, data.y, data.z, layersPermeabilities, x_dimensionLengths, y_dimensionLengths, heights);
            cartesianGrid.transmissibilityMultipliers = transmissibilitiyMultipliers;
            BaseBlock[] grid = cartesianGrid.Initialize();

            data.grid = grid;
        }

        private static void InitializePROPS(List<string> section, PVT pvt, SCAL scal, PorosityCalculator porosityCalculator)
        {

            var PROPS_KEYWORDS = new string[] { "SWFN", "SGFN", "SOF3", "ROCK", "DENSITY", "SO", "SG", "SW", "PRESSURE", "BUBBLEPOINT", "SWC", "PVTO", "PVTW", "PVTG", "SCAL"};

            double bubble_point;
            double[][] oil, oil_us, gas;
            oil = new double[4][]; oil_us = new double[4][]; gas = new double[4][];


            // rock reference pressure and compressibility
            double[] rock = Helper.GetData("ROCK", section, 2);
            porosityCalculator.pressure_ref = rock[0];
            porosityCalculator.Cf = rock[1];

            double[] densities = Helper.GetData("DENSITY", section, 3);

            So = Helper.GetData("SO", section, 1)[0];
            Sg = Helper.GetData("SG", section, 1)[0];
            Sw = Helper.GetData("SW", section, 1)[0];
            P = Helper.GetData("PRESSURE", section, 1)[0];
            Swc = Helper.GetData("SWC", section, 1)[0];

            bubble_point = Helper.GetData("BUBBLEPOINT", section, 1)[0];

            int index_pvto = section.FindIndex(x => x.Contains("PVTO"));
            int index_pvtg = section.FindIndex(x => x.Contains("PVTG"));
            int index_pvdg = section.FindIndex(x => x.Contains("PVDG"));
            int index_scal = section.FindIndex(x => x.Contains("SCAL"));

            int end_index = section.FindIndex(index_pvto + 1, x => x.ContainsAnyOf(PROPS_KEYWORDS));
            end_index = end_index < 0 ? section.Count : end_index;
            GetPVTO(oil, oil_us, section,index_pvto + 1, end_index);

            double[] water = Helper.GetData("PVTW", section, 5);

            if (index_pvdg >= 0)
            {
                end_index = section.FindIndex(index_pvdg + 1, x => x.ContainsAnyOf(PROPS_KEYWORDS));
                end_index = end_index < 0 ? section.Count : end_index;
                GetPVDG(gas, section, index_pvdg + 1, end_index);
            }

            if (index_pvtg >= 0)
            {
                end_index = section.FindIndex(index_pvtg + 1, x => x.ContainsAnyOf(PROPS_KEYWORDS));
                end_index = end_index < 0 ? section.Count : end_index;
                GetPVTG(gas, section, index_pvtg + 1, end_index);
            }

            //if (data.vaporizedOilPresent)
            //{
            //}
            //else
            //{
            //}

            int index_sgfn = section.FindIndex(x => x.Contains("SGFN"));
            end_index = section.FindIndex(index_sgfn + 1, x => x.ContainsAnyOf(PROPS_KEYWORDS));
            end_index = end_index < 0 ? section.Count : end_index;
            double[][] sgfn = new double[3][];
            GetSGFN(sgfn, section, index_sgfn + 1, end_index);

            int index_swfn = section.FindIndex(x => x.Contains("SWFN"));
            end_index = section.FindIndex(index_swfn + 1, x => x.ContainsAnyOf(PROPS_KEYWORDS));
            end_index = end_index < 0 ? section.Count : end_index;
            double[][] swfn = new double[3][];
            GetSWFN(swfn, section, index_swfn + 1, end_index);

            int index_sof3 = section.FindIndex(x => x.Contains("SOF3"));
            end_index = section.FindIndex(index_sof3 + 1, x => x.ContainsAnyOf(PROPS_KEYWORDS));
            end_index = end_index < 0 ? section.Count : end_index;
            double[][] sof3 = new double[3][];
            GetSOF3(sof3, section, index_sof3 + 1, end_index);


            pvt.Initialize(oil, oil_us, water, gas, bubble_point, densities, Swc);
            scal.Initialize(sgfn, swfn, sof3);

        }

        private static void GetPVDG(double[][] gas, List<string> section, int index_start, int index_end)
        {
            int count = index_end - index_start;
            gas[0] = new double[count];
            gas[1] = new double[count];
            gas[2] = new double[count];
            gas[3] = new double[count];

            for (int i = 0; i < count; i++)
            {
                var temp = Helper.GetData(section[i + index_start]);

                gas[0][i] = temp[0]; // pressure
                gas[3][i] = 0; // Rv
                gas[1][i] = temp[1]; // FVF
                gas[2][i] = temp[2]; // Viscosity
            }
        }

        private static void GetSOF3(double[][] sof3, List<string> section, int index_start, int index_end)
        {
            int count = index_end - index_start;
            sof3[0] = new double[count];
            sof3[1] = new double[count];
            sof3[2] = new double[count];

            for (int i = 0; i < count; i++)
            {
                var temp = Helper.GetData(section[i + index_start]);

                sof3[0][i] = temp[0]; // So
                sof3[1][i] = temp[1]; // Krow
                sof3[2][i] = temp[2]; // Krog
            }
        }

        private static void GetSWFN(double[][] swfn, List<string> section, int index_start, int index_end)
        {
            int count = index_end - index_start;
            swfn[0] = new double[count];
            swfn[1] = new double[count];
            swfn[2] = new double[count];

            for (int i = 0; i < count; i++)
            {
                var temp = Helper.GetData(section[i + index_start]);

                swfn[0][i] = temp[0]; // Sw
                swfn[1][i] = temp[1]; // Krw
                swfn[2][i] = temp[2]; // PCOW
            }
        }

        private static void GetSGFN(double[][] sgfn, List<string> section, int index_start, int index_end)
        {
            int count = index_end - index_start;
            sgfn[0] = new double[count];
            sgfn[1] = new double[count];
            sgfn[2] = new double[count];

            for (int i = 0; i < count; i++)
            {
                var temp = Helper.GetData(section[i + index_start]);

                sgfn[0][i] = temp[0]; // Sg
                sgfn[1][i] = temp[1]; // Krg
                sgfn[2][i] = temp[2]; // PCOG
            }
        }

        private static void GetPVTG(double[][] gas, List<string> section, int index_start, int index_end)
        {
            int count = index_end - index_start;
            gas[0] = new double[count];
            gas[1] = new double[count];
            gas[2] = new double[count];
            gas[3] = new double[count];

            for (int i = 0; i < count; i++)
            {
                var temp = Helper.GetData(section[i + index_start]);

                gas[0][i] = temp[0]; // pressure
                gas[3][i] = temp[1]; // Rv
                gas[1][i] = temp[2]; // FVF
                gas[2][i] = temp[3]; // Viscosity
            }
        }

        private static void GetPVTO(double[][] oil, double[][] oil_us, List<string> section, int start_index, int end_index)
        {
            int highest_rs_index = 0;
            int counter = 0;

            var temp_array = new double[end_index - start_index][];
            for (int i = 0; i < temp_array.Length; i++)
            {
                temp_array[i] = Helper.GetData(section[i + start_index]);
            }

            int saturated_count = temp_array.Count(x => x.Length == 4);
            int under_saturated_count = temp_array.Length - saturated_count + 1;

            oil[0] = new double[saturated_count];
            oil[1] = new double[saturated_count];
            oil[2] = new double[saturated_count];
            oil[3] = new double[saturated_count];

            oil_us[0] = new double[under_saturated_count];
            oil_us[1] = new double[under_saturated_count];
            oil_us[2] = new double[under_saturated_count];
            oil_us[3] = new double[under_saturated_count];

            for (int i = 0; i < temp_array.Length; i++)
            {
                var temp = temp_array[i];

                if (temp.Length == 3)
                {
                    if (highest_rs_index == 0)
                    {
                        highest_rs_index = counter - 1;
                        counter = 0;

                        oil_us[0][counter] = oil[0][highest_rs_index];
                        oil_us[1][counter] = oil[1][highest_rs_index];
                        oil_us[2][counter] = oil[2][highest_rs_index];
                        oil_us[3][counter] = oil[3][highest_rs_index];
                        counter++;
                    }

                    oil_us[0][counter] = temp[0];
                    oil_us[1][counter] = temp[1];
                    oil_us[2][counter] = temp[2];
                    oil_us[3][counter] = oil[3][highest_rs_index];
                }
                else // 4 columns
                {
                    oil[0][counter] = temp[1];
                    oil[1][counter] = temp[2];
                    oil[2][counter] = temp[3];
                    oil[3][counter] = temp[0];
                }

                counter++;
            }
        }

        private static void InitializeSOLUTION(List<string> section)
        {
            // To-Do
        }

        private static void InitializeSUMMARY(List<string> section)
        {
            string[] wellKeyWords = new string[] { "WOPR", "WWPR", "WGPR", "WGPRS", "WGPRF", "WBHP", "WGOR" };
            string[] blockKeyWords = new string[] { "BOFVF", "BWFVF", "BGFVF", "BOSAT", "BWSAT", "BGSAT", "BKRO", "BKRW", "BKRG", "BPR", "BRS" };
            string[] singleKeyWords = new string[] { "CONS", "FILE", "FOPR", "MBEO", "MBEW", "MBEG", "NEWTON", "TCPUTS", "FGIP", "FGIPL", "FGIPG"};

            string[] allKeyWords = new string[wellKeyWords.Length + blockKeyWords.Length + singleKeyWords.Length];
            Array.Copy(wellKeyWords, allKeyWords, wellKeyWords.Length);
            Array.Copy(blockKeyWords, 0, allKeyWords, wellKeyWords.Length, blockKeyWords.Length);
            Array.Copy(singleKeyWords, 0, allKeyWords, wellKeyWords.Length + blockKeyWords.Length, singleKeyWords.Length);

            int index;
            string[] wells;
            int[] indices;
            string[] xyzIndices;
            //string keyword;

            output.console = section.FindIndex(line => line.Contains("CONS")) != -1 ? true : false;
            output.file = section.FindIndex(line => line.Contains("FILE")) != -1 ? true : false;

            foreach (string keyword in wellKeyWords)
            {
                index = section.FindIndex(line => line.Contains(keyword));
                // make sure the keyword exists
                if (index != -1)
                {
                    var temp = new List<string>();
                    while (index < section.Count - 1 && !wellKeyWords.Contains(section[index + 1]))
                    {
                        temp.AddRange(Helper.GetDataAsString(section[index + 1]));
                        index++;
                    }
                    wells = temp.ToArray();
                    indices = data.wells.Where(x => x.name.ContainsAnyOf(wells)).Select(x => x.index).ToArray();
                    output.wellKeyWords.Add(keyword);
                    output.wellsIndices.Add(indices);
                }

            }


            foreach (string keyword in blockKeyWords)
            {
                index = section.FindIndex(line => line.Contains(keyword));

                if (index != -1)
                {
                    int[][] indices_array = Helper.GetMultipleRowsIndices(keyword, allKeyWords, section);
                    indices = new int[indices_array.Length];
                    xyzIndices = new string[indices_array.Length];

                    for (int i = 0; i < indices_array.Length; i++)
                    {
                        indices[i] = Misc.Rectangular.xyzToNatural(data.x, data.y, data.z, indices_array[i][0], indices_array[i][1], indices_array[i][2]);
                    }

                    for (int i = 0; i < indices_array.Length; i++)
                    {
                        xyzIndices[i] = $"{indices_array[i][0] + 1},{indices_array[i][1] + 1},{indices_array[i][2] + 1}";
                    }

                    output.blockKeyWords.Add(keyword);
                    output.blocksIndices.Add(indices);
                    output.blocksXYZIndices.Add(xyzIndices);
                }

            }

            for (int i = 2; i < singleKeyWords.Length; i++)
            {
                index = section.FindIndex(line => line.Contains(singleKeyWords[i]));

                if (index != -1)
                {
                    output.singleKeyWords.Add(singleKeyWords[i]);
                }
            }

        }

        private static void InitializeSCHEDULE(SimulationData data, List<string> section)
        {
            int x, y, z, well_index;
            double well_radius, skin_factor;
            Global.WellType type;
            Global.WellControl control = Global.WellControl.OilRate;
            double flow_rate, BHP;

            int index = section.FindIndex(line => line.Contains("WELSPECS"));
            string[][] WELSPECS = new string[data.wells.Length][];

            for (int i = 0; i < data.wells.Length; i++)
            {
                WELSPECS[i] = Helper.GetDataAsString(section[index + 1 + i]);
            }

            

            index = section.FindIndex(line => line.Contains("WELCONT"));
            string[][] WELCONT = new string[data.wells.Length][];

            for (int i = 0; i < data.wells.Length; i++)
            {
                WELCONT[i] = Helper.GetDataAsString(section[index + 1 + i]);
            }


            for (int i = 0; i < data.wells.Length; i++)
            {

                x = int.Parse(WELSPECS[i][1]) - 1;
                y = int.Parse(WELSPECS[i][2]) - 1;
                z = int.Parse(WELSPECS[i][3]) - 1;

                well_index = Misc.Rectangular.xyzToNatural(data.x, data.y, data.z, x, y, z);

                skin_factor = double.Parse(WELSPECS[i][4]);

                well_radius = double.Parse(WELSPECS[i][5]);

                //type = WELCONT[i][0].Contains("PRODUCER") ? Global.WellType.Production : Global.WellType.Injection;
                int well_control_index;
                for (well_control_index = 0; well_control_index < WELCONT.Length; well_control_index++)
                {
                    if (WELCONT[well_control_index][0] == WELSPECS[i][0]/*well name*/)
                    {
                        break;
                    }
                }
                type = WELCONT[well_control_index][4].Contains("PRODUCER") ? Global.WellType.Production : Global.WellType.Injection;

                if (type == Global.WellType.Production)
                {
                    control = WELCONT[well_control_index][1] == "RATE" ? Global.WellControl.OilRate : Global.WellControl.BHP;
                }
                else if (type == Global.WellType.Injection)
                {
                    // RATE or GRATE means gas injection
                    // WRATE means water injection
                    string temp = WELCONT[well_control_index][1];
                    control = temp == "RATE" || temp == "GRATE" ? Global.WellControl.GasRate : temp == "WRATE" ? Global.WellControl.WaterRate : Global.WellControl.BHP;
                }

                flow_rate = double.Parse(WELCONT[well_control_index][2]);
                BHP = double.Parse(WELCONT[well_control_index][3]);

                data.wells[i] = new Well.BaseWell(data, well_index, type, control, well_radius, skin_factor, BHP, flow_rate);

                data.wells[i].name = WELSPECS[i][0];

            }

            data.originalTimeStep = Helper.GetData("TSTEP", section)[0];
            data.endingTime = Helper.GetData("ENDTIME", section)[0];
        }

        // get data from initialization file.

        private static void InitializationFile(List<string> lines)
        {
            data.timeStepSlashingFactor = Helper.GetData("TSTEPFACTOR", lines, 1)[0];
            data.originalRelaxationFactor = Helper.GetData("RELAXFACTOR", lines, 1)[0];
            data.minimumRelaxation = Helper.GetData("MINRELAX", lines, 1)[0];
            data.relaxationFactorDecrement = Helper.GetData("RELAXDEC", lines, 1)[0];
            data.maximumNonLinearIterations = Helper.GetData("MAXNONLINIT", lines, 1)[0];
            data.maximumMaterialBalanceErrorRatio = Helper.GetData("MBERATIO", lines, 1)[0];
            data.MBE_Tolerance = Helper.GetData("MAXMBE", lines, 1)[0];
            Global.padding = (int)Helper.GetData("PADDING", lines, 1)[0];
            Global.decimalPlaces = "#0." + new string('0', (int)Helper.GetData("DECIPLCS", lines, 1)[0]);
            Global.EPSILON_P = Helper.GetData("EPSILON_P", lines, 1)[0];
            Global.EPSILON_S = Helper.GetData("EPSILON_S", lines, 1)[0];
            Global.MINIMUM = Helper.GetData("DERIVATIVE_MINIMUM", lines, 1)[0];
        }


        // helper method for adding each input data file section to its corresponding list.
        private static void AddSection(string keyWord, List<string> section, List<string> lines, string[] dataFileSections)
        {

            for (int i = endIndex; i < lines.Count; i++)
            {
                if (!startIndexSet && lines[i].Contains(keyWord))
                {
                    startIndex = i;
                    startIndexSet = true;
                }

                if (i != startIndex && lines[i].ContainsAnyOf(dataFileSections))
                {
                    endIndex = i;
                    startIndexSet = false;
                    break;
                }
            }

            for (int i = startIndex + 1; i < endIndex; i++)
            {
                section.Add(lines[i]);
            }
        }

    }
}