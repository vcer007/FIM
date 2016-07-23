/// <summary>
/// This namespace contains helpful miscellaneous methods that are not re;ated to an essential ctegory.
/// </summary>
namespace FIM.Misc
{
    /// <summary>
    /// This class contains helper method for dealing with 3D cartesian grids.
    /// </summary>
    /// <remarks>
    /// As FIM deals with PeBi grids, it is customary to convert any type of grids to a linear array of blocks.
    /// This is possible by assigning <see cref="FIM.Core.BaseBlock.neighbourBlocksIndices"/> to each <see cref="FIM.Core.BaseBlock"/>.
    /// </remarks>
    abstract class Rectangular
    {
        /// <summary>
        /// Convert cartesian 3D dimensions to natural "liear" dimensions.
        /// </summary>
        /// <param name="x">The x.</param>
        /// <param name="y">The y.</param>
        /// <param name="z">The z.</param>
        /// <param name="i">The i.</param>
        /// <param name="j">The j.</param>
        /// <param name="k">The k.</param>
        /// <returns>The index in the grid</returns>
        /// <seealso cref="FIM.Core.SimulationData.grid"/>
        public static int xyzToNatural(int x, int y, int z, int i, int j, int k)
        {
            return i + j * x + k * x * y;
        }
    }
}
