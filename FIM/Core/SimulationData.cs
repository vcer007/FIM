using FIM.Fluid;
using FIM.Rock;
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
        public double minimum_relaxation = 0.5;

        public double original_relaxation_factor = 1;
        public double original_time_step;

        public double relaxation_factor = 1;
        public double time_step;


        public double MBE_Gas;


        public int x, y, z;

        public Global.Phase[] phases;
        public bool solubleGasPresent;

        public double tolerance;

        // PVT and Rock data calculators
        public PVT pvt;
        public Kr kr;
        public Porosity porosity;

        // the grid "a list of the blocks"
        public BaseBlock[] grid;

        public SimulationData(int x, int y, int z, BaseBlock[] grid)
        {
            this.x = x; this.y = y; this.z = z;

            this.grid = grid;
        }
    }
}
