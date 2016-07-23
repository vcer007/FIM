using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FIM.FluidData;
using FIM.Core;
using FIM.Solver;

namespace FIM
{
    class Program
    {
        static void Main(string[] args)
        {
            SimulationData data = Initialization.Model.initiaize();
            FullyImplicitSolver.RunSimulation(data);

            Console.ReadKey();
        }
    }
}