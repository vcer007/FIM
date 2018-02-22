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

        static double So, Sw, Sg, P;

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

            #region Initialization File

            lines = Helper.ReadFile(initializationFilePath);

            InitializationFile(lines);

            //data.timeStepSlashingFactor = 0.5;
            //data.originalRelaxationFactor = 1;
            //data.minimumRelaxation = 0.5;
            //data.relaxationFactorIncrement = -0.1;
            //data.maximumNonLinearIterations = 25;
            //data.maximumConvergenceErrorRatio = 0.5;
            //data.MBE_Tolerance = 1;

            //data.relaxationFactor = data.originalRelaxationFactor;
            //data.endingTime = 10 * 365;
            //data.originalTimeStep = 1;

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
                //phases.Add(Global.Phase.DissolvedGas);
                data.solubleGasPresent = true;
            }

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
            double bubble_point;
            double[][] oil, oil_us, water, water_us, gas;
            oil = new double[5][]; oil_us = new double[5][]; water = new double[5][]; water_us = new double[5][]; gas = new double[4][];
            double[][] Kr_data = new double[4][];


            // rock reference pressure and compressibility
            double[] rock = Helper.GetData("ROCK", section, 2);
            porosityCalculator.pressure_ref = rock[0];
            porosityCalculator.Cf = rock[1];

            So = Helper.GetData("SO", section, 1)[0];
            Sg = Helper.GetData("SG", section, 1)[0];
            Sw = Helper.GetData("SW", section, 1)[0];
            P = Helper.GetData("PRESSURE", section, 1)[0];

            bubble_point = Helper.GetData("BUBBLEPOINT", section, 1)[0];

            int index = section.FindIndex(x => x.Contains("PVTO"));
            oil[0] = Helper.GetData(section[index + 1]);
            oil[1] = Helper.GetData(section[index + 2]);
            oil[2] = Helper.GetData(section[index + 3]);
            oil[3] = Helper.GetData(section[index + 4]);
            oil[4] = Helper.GetData(section[index + 5]);

            oil_us[0] = Helper.GetData(section[index + 6]);
            oil_us[1] = Helper.GetData(section[index + 7]);
            oil_us[2] = Helper.GetData(section[index + 8]);
            oil_us[3] = Helper.GetData(section[index + 9]);
            oil_us[4] = Helper.GetData(section[index + 10]);

            index = section.FindIndex(x => x.Contains("PVTW"));

            water[0] = Helper.GetData(section[index + 1]);
            water[1] = Helper.GetData(section[index + 2]);
            water[2] = Helper.GetData(section[index + 3]);
            water[3] = Helper.GetData(section[index + 4]);
            water[4] = Helper.GetData(section[index + 5]);

            water_us[0] = Helper.GetData(section[index + 6]);
            water_us[1] = Helper.GetData(section[index + 7]);
            water_us[2] = Helper.GetData(section[index + 8]);
            water_us[3] = Helper.GetData(section[index + 9]);
            water_us[4] = Helper.GetData(section[index + 10]);


            index = section.FindIndex(x => x.Contains("PVDG"));

            gas[0] = Helper.GetData(section[index + 1]);
            gas[1] = Helper.GetData(section[index + 2]);
            gas[2] = Helper.GetData(section[index + 3]);
            gas[3] = Helper.GetData(section[index + 4]);


            index = section.FindIndex(x => x.Contains("SCAL"));

            Kr_data[0] = Helper.GetData(section[index + 1]);
            Kr_data[1] = Helper.GetData(section[index + 2]);
            Kr_data[2] = Helper.GetData(section[index + 3]);
            Kr_data[3] = Helper.GetData(section[index + 4]);

            pvt.Initialize(oil, oil_us, water, water_us, gas, bubble_point);
            scal.Initialize(Kr_data);

        }

        private static void InitializeSOLUTION(List<string> section)
        {
            // To-Do
        }

        private static void InitializeSUMMARY(List<string> section)
        {
            string[] wellKeyWords = new string[] { "WOPR", "WWPR", "WGPR", "WGPRS", "WGPRF", "WBHP", "WGOR" };
            string[] blockKeyWords = new string[] { "BOFVF", "BWFVF", "BGFVF", "BOSAT", "BWSAT", "BGSAT", "BKRO", "BKRW", "BKRG", "BPR", "BRS" };
            string[] singleKeyWords = new string[] { "CONS", "FILE", "FOPR", "MBEO", "MBEW", "MBEG", "NEWTON", "TCPUTS"};

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
                    control = WELCONT[well_control_index][1] == "RATE" ? Global.WellControl.GasRate : Global.WellControl.BHP;
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