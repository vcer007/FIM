using FIM.Core;
using FIM.FluidData;
using FIM.Initialization;
using FIM.RockData;
using FIM.Extensions;

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace FIM.Parser
{
    class Eclipse
    {
        static string[] dataFileSections = new string[] { "RUNSPEC", "GRID", "PROPS", "SOLUTION", "SUMMARY", "SCHEDULE" , "END"};

        static List<string> RUNSPEC_section = new List<string>();
        static List<string> GRID_section = new List<string>();
        static List<string> PROPS_section = new List<string>();
        static List<string> SOLUTION_section = new List<string>();
        static List<string> SUMMARY_section = new List<string>();
        static List<string> SCHEDULE_section = new List<string>();

        static int startIndex = 0, endIndex = 0;
        static bool startIndexSet;

        static SimulationData data = new SimulationData();
        static PorosityCalculator porosityCalculator = new PorosityCalculator();
        static PVT pvt = new PVT();
        static SCAL scal = new SCAL();

        static double So, Sw, Sg, P;

        public static SimulationData ReadInputData(string dataFilePath, string initializationFilePath)
        {
            List<string> lines;


            #region Data File

            lines = ReadFile(dataFilePath);

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
            InitializeSUMMARY(SUMMARY_section);
            InitializeSCHEDULE(data, SCHEDULE_section);

            data.porosityCalculator = porosityCalculator;
            data.pvt = pvt;
            data.scal = scal;


            // update block properties
            for (int i = 0; i < data.grid.Length; i++)
            {
                data.grid[i].UpdateProperties(data, P, Sw, Sg, So, 0);
            }

            #endregion

            #region Initialization File

            lines = ReadFile(initializationFilePath);

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

        private static void InitializationFile(List<string> lines)
        {
            data.timeStepSlashingFactor = getData("TSTEPFACTOR", lines, 1)[0];
            data.originalRelaxationFactor = getData("RELAXFACTOR", lines, 1)[0];
            data.minimumRelaxation = getData("MINRELAX", lines, 1)[0];
            data.relaxationFactorDecrement = getData("RELAXDEC", lines, 1)[0];
            data.maximumNonLinearIterations = getData("MAXNONLINIT", lines, 1)[0];
            data.maximumConvergenceErrorRatio = getData("MAXCONV", lines, 1)[0];
            data.MBE_Tolerance = getData("MAXMBE", lines, 1)[0];
        }

        private static List<string> ReadFile(string filePath)
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
            x_dimensionLengths = getData(section[index + 1], data.x);
            // y
            index = section.FindIndex(x => x.Contains("DY"));
            y_dimensionLengths = getData(section[index + 1], data.y);
            // z
            index = section.FindIndex(x => x.Contains("DZ"));
            heights = getData(section[index + 1], data.z);

            // permeabilities

            // permeability in the x direction
            index = section.FindIndex(x => x.Contains("PERMX"));
            layersPermeabilities = getData(section[index + 1], data.z);

            // transmissibilities multipliers in the z direction only.
            index = section.FindIndex(x => x.Contains("MULTZ"));
            if (index != -1)
            {
                transmissibilitiyMultipliers = getData(section[index + 1], data.z);
            }
            else
            {
                transmissibilitiyMultipliers = Enumerable.Repeat(1.0, data.z).ToArray();
            }

            // porosity
            index = section.FindIndex(x => x.Contains("PORO"));
            // only one value of porosity is initially assigned to the whole model.
            porosityCalculator.porosity_ref = getData(section[index + 1], 1)[0];



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
            double[] rock = getData("ROCK", section, 2);
            porosityCalculator.pressure_ref = rock[0];
            porosityCalculator.Cf = rock[1];

            So = getData("SO", section, 1)[0];
            Sg = getData("SG", section, 1)[0];
            Sw = getData("SW", section, 1)[0];
            P = getData("PRESSURE", section, 1)[0];

            bubble_point = getData("BUBBLEPOINT", section, 1)[0];

            int index = section.FindIndex(x => x.Contains("PVTO"));
            oil[0] = getData(section[index + 1]);
            oil[1] = getData(section[index + 2]);
            oil[2] = getData(section[index + 3]);
            oil[3] = getData(section[index + 4]);
            oil[4] = getData(section[index + 5]);

            oil_us[0] = getData(section[index + 6]);
            oil_us[1] = getData(section[index + 7]);
            oil_us[2] = getData(section[index + 8]);
            oil_us[3] = getData(section[index + 9]);
            oil_us[4] = getData(section[index + 10]);

            index = section.FindIndex(x => x.Contains("PVTW"));

            water[0] = getData(section[index + 1]);
            water[1] = getData(section[index + 2]);
            water[2] = getData(section[index + 3]);
            water[3] = getData(section[index + 4]);
            water[4] = getData(section[index + 5]);

            water_us[0] = getData(section[index + 6]);
            water_us[1] = getData(section[index + 7]);
            water_us[2] = getData(section[index + 8]);
            water_us[3] = getData(section[index + 9]);
            water_us[4] = getData(section[index + 10]);


            index = section.FindIndex(x => x.Contains("PVDG"));

            gas[0] = getData(section[index + 1]);
            gas[1] = getData(section[index + 2]);
            gas[2] = getData(section[index + 3]);
            gas[3] = getData(section[index + 4]);


            index = section.FindIndex(x => x.Contains("SCAL"));

            Kr_data[0] = getData(section[index + 1]);
            Kr_data[1] = getData(section[index + 2]);
            Kr_data[2] = getData(section[index + 3]);
            Kr_data[3] = getData(section[index + 4]);

            pvt.Initialize(oil, oil_us, water, water_us, gas, bubble_point);
            scal.Initialize(Kr_data);

        }

        private static void InitializeSOLUTION(List<string> section)
        {
            // To-Do
        }

        private static void InitializeSUMMARY(List<string> section)
        {
            // To-Do
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
                WELSPECS[i] = getDataAsString(section[index + 1 + i]);
            }

            

            index = section.FindIndex(line => line.Contains("WELCONT"));
            string[][] WELCONT = new string[data.wells.Length][];

            for (int i = 0; i < data.wells.Length; i++)
            {
                WELCONT[i] = getDataAsString(section[index + 1 + i]);
            }


            for (int i = 0; i < data.wells.Length; i++)
            {

                x = int.Parse(WELSPECS[i][1]) - 1;
                y = int.Parse(WELSPECS[i][2]) - 1;
                z = int.Parse(WELSPECS[i][3]) - 1;

                well_index = Misc.Rectangular.xyzToNatural(data.x, data.y, data.z, x, y, z);

                skin_factor = double.Parse(WELSPECS[i][4]);

                well_radius = double.Parse(WELSPECS[i][5]);

                type = WELCONT[i][0].Contains("PRODUCER") ? Global.WellType.Production : Global.WellType.Injection;

                if (type == Global.WellType.Production)
                {
                    control = WELCONT[i][1] == "RATE" ? Global.WellControl.OilRate : Global.WellControl.BHP;
                }
                else if (type == Global.WellType.Injection)
                {
                    control = WELCONT[i][1] == "RATE" ? Global.WellControl.GasRate : Global.WellControl.BHP;
                }

                flow_rate = double.Parse(WELCONT[i][2]);
                BHP = double.Parse(WELCONT[i][3]);

                data.wells[i] = new Well.BaseWell(data, well_index, type, control, well_radius, skin_factor, BHP, flow_rate);

                data.wells[i].name = WELSPECS[i][0];

            }

            data.originalTimeStep = getData("TSTEP", section)[0];
            data.endingTime = getData("ENDTIME", section)[0];
        }




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

        private static double[] getData(string line, int size = 1)
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

        private static double[] getData(string keyWord, List<string> section, int size = 1)
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

        private static string[] getDataAsString(string line, int size = 1)
        {
            string[] data = Regex.Split(line, @"\b\s+\b");

            return data;
        }

    }
}