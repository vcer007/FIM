using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FIM.Fluid;
using FIM.Core;
using FIM.Solver;

namespace FIM
{
    class Program
    {
        static void Main(string[] args)
        {
            SimulationData data = Initialize.Odeh.initiaize();
            double[] R = FullyImplicit.calculate_R(data);
            double[][] jacobians = FullyImplicit.calculateJacobians(data, R);

            Console.WriteLine("Done");
        }
    }
}
