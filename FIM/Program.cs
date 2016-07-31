using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FIM.FluidData;
using FIM.Core;
using FIM.Solver;
using FIM.Initialization;

namespace FIM
{
    /// <summary>
    /// The starting point of the simulator
    /// </summary>
    public class Program
    {
        /// <summary>
        /// The first method executed in the simulator.
        /// </summary>
        /// <remarks>
        /// C# runs this method automatically at start.
        /// </remarks>
        public static void Main()
        {
            //SimulationData data = Initialization.Model.initiaize();

            // display a instructive intro.
            WriteIntro();

            // get the initialization file name as an input.
            Console.WriteLine("Enter the initialization file name : ");
            string initialization_file_name = Console.ReadLine();

            // get the input data file name as an input.
            Console.WriteLine("Enter the input data file name : ");
            string data_file_name = Console.ReadLine();

            Console.WriteLine();

            // read the initialization file and the input data file.
            SimulationData data = Parser.Eclipse.ReadInputData(data_file_name, initialization_file_name);

            // display some info read from the input data file.
            WriteInfo(data);

            // determine which solver is to be used.
            if (data.solutionProcedure == Global.SolutionProcedure.FullyImplicit)
            {
                FullyImplicitSolver.RunSimulation(data);
            }
            else
            {
                IMPES_Solver.RunSimulation(data);
            }


            // this prevents the console window from automatically closing after 
            Console.ReadKey();
        }

        private static void WriteIntro()
        {
            Console.WriteLine(".................FIM Black-Oil Simulator..................");
            Console.WriteLine();
            Console.WriteLine("1- Make Sure the files are located within the same directory as the simulator.");
            Console.WriteLine("2- Spell the file names correctly with their respective extensions, '.INIT' and '.DATA'.");
            Console.WriteLine("3- Names are case insensitive.");
            Console.WriteLine();
        }

        private static void WriteInfo(SimulationData data)
        {
            Console.WriteLine();

            Console.Write("Title         : ");
            Console.Write(data.title);
            Console.WriteLine();

            Console.Write("Solver        : ");
            if (data.solutionProcedure == Global.SolutionProcedure.FullyImplicit)
            {
                Console.Write("Fully Implicit");
            }
            else
            {
                Console.Write("IMPES");
            }
            Console.WriteLine();

            Console.Write("Ending Time   : ");
            Console.Write(data.endingTime + " ");
            Console.Write("days");
            Console.WriteLine();

            Console.Write("Grid Size     : ");
            Console.Write(data.grid.Length + " ");
            Console.Write("blocks");
            Console.WriteLine();

            Console.Write("Wells Present : ");
            Console.Write(data.wells.Length + " ");
            Console.Write("wells");
            Console.WriteLine();

            Console.WriteLine();
            Console.WriteLine();

            Console.WriteLine("Press any key to begin the simulation ...");
            Console.ReadKey();

            Console.WriteLine();
        }
    }
}