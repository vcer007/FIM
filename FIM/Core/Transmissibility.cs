using System;

namespace FIM.Core
{
    /// <summary>
    /// This class is used to calculate the geometrical transmissibility between two blocks.
    /// </summary>
    public static class Transmissibility
    {
        /// <summary>
        /// Calculates the geometrical transmissibility between two blocks.
        /// </summary>
        /// <remarks>
        /// This method is based on harmonic averaging.
        /// The geometrical transmissibility factor calculated from this method is fixed throughout the simulation.
        /// It's independent of pressures and saturations.
        /// </remarks>
        /// <param name="block_1">block_1.</param>
        /// <param name="block_2">block_2.</param>
        /// <param name="multiplier">a multiplier used for altering inter-layers transmissibilities.</param>
        /// <returns>the Geometrical transmissibility between two blocks.</returns>
        /// <seealso cref="BaseBlock"/>
        public static double Calculate(BaseBlock block_1, BaseBlock block_2, double multiplier = 1)
        {
            // As each block stores a list of side faces areas, we need to find the index of each block relative to the other.
            int index_1, index_2;

            index_1 = Array.IndexOf(block_1.neighborBlocksIndices, block_2.index);
            index_2 = Array.IndexOf(block_2.neighborBlocksIndices, block_1.index);

            double length_1, area_1, permeability_1;
            double length_2, area_2, permeability_2;

            length_1 = block_1.deltaXList[index_1];
            area_1 = block_1.areaList[index_1];
            permeability_1 = block_1.permeability[index_1];

            length_2 = block_2.deltaXList[index_2];
            area_2 = block_2.areaList[index_2];
            permeability_2 = block_2.permeability[index_2];

            double G = Global.Bc / (length_1 / (area_1 * permeability_1) + length_2 / (area_2 * permeability_2));
            return G * multiplier;
        }
    }
}
