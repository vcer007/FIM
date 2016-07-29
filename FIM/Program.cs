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

            SimulationData data = Parser.Eclipse.ReadFile("Odeh.Data");

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
    }
}