using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Extensions.IMPES
{
    public static class Compressibility
    {
        public static double GetCompressibilityTotal(this BaseBlock block, SimulationData data)
        {
            return GetCompressibilityOil(block) * block.So[0] + GetCompressibilityWater(block) * block.Sw[0] + GetCompressibilityGas(block) * block.Sg[0] + data.porosityCalculator.Cf;
        }

        private static double GetCompressibilityOil(BaseBlock block)
        {
            return -1 / block.Bo[0] * (block.Bo[1] - block.Bo[0]) / (block.P[1] - block.P[0]) + block.Bg[1] / block.Bo[0] * (block.Rso[1] - block.Rso[0]) / (block.P[1] - block.P[0]);
        }

        private static double GetCompressibilityGas(BaseBlock block)
        {
            return -1 / block.Bg[0] * (block.Bg[1] - block.Bg[0]) / (block.P[1] - block.P[0]);
        }

        private static double GetCompressibilityWater(BaseBlock block)
        {
            return -1 / block.Bw[0] * (block.Bw[1] - block.Bw[0]) / (block.P[1] - block.P[0]);
        }

        
    }
}
