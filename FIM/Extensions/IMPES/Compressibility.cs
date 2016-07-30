using FIM.Core;

/// <summary>
/// This namespace contains extension methods used for the IMPES solver.
/// </summary>
namespace FIM.Extensions.IMPES
{
    /// <summary>
    /// This class contains methods to calculate fluids compressibilities.
    /// </summary>
    public static class Compressibility
    {
        /// <summary>
        /// Gets the total compressibility.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <param name="data">The <see cref="SimulationData"/> data.</param>
        /// <returns></returns>
        public static double GetCompressibilityTotal(this BaseBlock block, SimulationData data)
        {
            return GetCompressibilityOil(block) * block.So[0] + GetCompressibilityWater(block) * block.Sw[0] + GetCompressibilityGas(block) * block.Sg[0] + data.porosityCalculator.Cf;
        }

        /// <summary>
        /// Gets the compressibility of oil.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <returns></returns>
        private static double GetCompressibilityOil(BaseBlock block)
        {
            return -1 / block.Bo[0] * (block.Bo[1] - block.Bo[0]) / (block.P[1] - block.P[0]) + block.Bg[1] / block.Bo[0] * (block.Rso[1] - block.Rso[0]) / (block.P[1] - block.P[0]);
        }

        /// <summary>
        /// Gets the compressibility of gas.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <returns></returns>
        private static double GetCompressibilityGas(BaseBlock block)
        {
            return -1 / block.Bg[0] * (block.Bg[1] - block.Bg[0]) / (block.P[1] - block.P[0]);
        }

        /// <summary>
        /// Gets the compressibility of water.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <returns></returns>
        private static double GetCompressibilityWater(BaseBlock block)
        {
            return -1 / block.Bw[0] * (block.Bw[1] - block.Bw[0]) / (block.P[1] - block.P[0]);
        }

        
    }
}
