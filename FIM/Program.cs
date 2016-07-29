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
    class Program
    {
        static void Main(string[] args)
        {
            //SimulationData data = Initialization.Model.initiaize();

            Console.WriteLine(".................FIM Black-Oil Simulator..................");
            Console.WriteLine("Make Sure the files are located within the same directory as the simulator .exe");
            Console.WriteLine("Spell the file names correctly with their respective extensions, .INIT and .DATA.");
            Console.WriteLine("Names are case insensitive.");
            Console.WriteLine();

            Console.WriteLine("Enter the initialization file name : ");
            string initialization_file_name = Console.ReadLine();

            Console.WriteLine("Enter the input data file name : ");
            string data_file_name = Console.ReadLine();

            Console.WriteLine();

            SimulationData data = Parser.Eclipse.ReadInputData(data_file_name, initialization_file_name);



            WriteInfo(data);



            if (data.solutionProcedure == Global.SolutionProcedure.FullyImplicit)
            {
                FullyImplicitSolver.RunSimulation(data);
            }
            else
            {
                IMPES_Solver.RunSimulation(data);
            }

            Console.ReadKey();
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
        }
    }
}