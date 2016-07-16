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
            FullyImplicit.RunSimulation(data);

            //SimulationData data = Initialize.ThreeBlocksLinear.initiaize();
            //FIM_ThreeBlocksLinear.RunSimulation(data);

            //SimulationData data = Initialize.SingleBlockModel.initiaize();
            //FIM_block.RunSimulation(data);

            //SimulationData data = Initialize.SingleLayerModel.initiaize();
            //FIM_SingleLayer.RunSimulation(data);

            Console.ReadKey();
        }
    }
}