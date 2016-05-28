using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Core
{
    /// <summary>
    /// This class contains all the general data of the simulation run "Run Specifications".
    /// </summary>
    class SimulationData
    {
        public double time_stpe;

        public int x, y, z;

        public Global.Phase[] phases;
        public bool solubleGasPresent;

        // the grid "a list of the blocks"
        public BaseBlock[] grid;

        public SimulationData(int x, int y, int z, BaseBlock[] grid)
        {
            this.x = x; this.y = y; this.z = z;

            this.grid = grid;
        }
    }
}
