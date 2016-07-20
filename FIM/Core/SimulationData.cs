using FIM.FluidData;
using FIM.RockData;
using FIM.Well;
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
    public class SimulationData
    {
        /// <summary>
        /// The minimum relaxation factor value.
        /// </summary>
        /// <remarks>
        /// <para>Default value is equal to 0.5.</para>
        /// <para>If the solution does not converge till this value, the time step is reduced.</para>
        /// </remarks>
        /// <seealso cref="SimulationData.time_step_slashing_factor"/>
        /// <seealso cref="SimulationData.relaxation_factor"/>
        public double minimum_relaxation = 0.5;

        /// <summary>
        /// The original relaxation factor value.
        /// </summary>
        /// <remarks>
        /// The <see cref="relaxation_factor"/> is set to this value at the beginning of the time step.
        /// </remarks>
        /// <seealso cref="SimulationData.minimum_relaxation"/>
        public double original_relaxation_factor = 1;

        /// <summary>
        /// The original_time_step
        /// </summary>
        public double original_time_step;

        /// <summary>
        /// The relaxation factor used when the solution oscillates.
        /// </summary>
        /// <remarks>
        /// When oscillations are detected, this value is used to "dampen" the calculated dx.
        /// </remarks>
        /// <seealso cref="original_relaxation_factor"/>
        /// <seealso cref="minimum_relaxation"/>
        public double relaxation_factor = 1;

        /// <summary>
        /// The time step used in calculations.
        /// </summary>
        public double time_step;

        /// <summary>
        /// The factor used to reduce the <see cref="SimulationData.time_step"/> when the solver stagnates or oscillates.
        /// </summary>
        public double time_step_slashing_factor = 0.5;

        /// <summary>
        /// The material balance error between fluid originally and currently in-place.
        /// </summary>
        /// <seealso cref="tolerance"/>
        public double MBE_Gas, MBE_Oil;

        /// <summary>
        /// The phases available in the model.
        /// </summary>
        /// <seealso cref="Global.Phase"/>
        public Global.Phase[] phases;

        /// <summary>
        /// A boolean to indicate if soluble gas is present.
        /// </summary>
        /// <seealso cref="phases"/>
        public bool solubleGasPresent;

        /// <summary>
        /// The tolerance used to determine convergence.
        /// </summary>
        /// <remarks>
        /// currently tolerance is set as an absolute variance between original fluid in place and fluid in place.
        /// This behavior is used to ensure very high accuracy.
        /// However, a percentage may be used to yield sufficiently accurate solutions with a lower number of iterations required.
        /// </remarks>
        public double tolerance;

        /// <summary>
        /// An instance of the <see cref="PVT"/> class.
        /// </summary>
        public PVT pvt;

        /// <summary>
        /// An instance of the <see cref="SCAL"/> class.
        /// </summary>
        public SCAL scal;

        /// <summary>
        /// An instance of the <see cref="PorosityCalculator"/> class.
        /// </summary>
        public PorosityCalculator porosity_calculator;

        // the grid "a list of the blocks"        
        /// <summary>
        /// An array of <see cref="BaseBlock"/> that compose the grid.
        /// </summary>
        public BaseBlock[] grid;

        public WellData[] wells;

        /// <summary>
        /// Initializes a new instance of the <see cref="SimulationData"/> class.
        /// </summary>
        /// <param name="grid">The <see cref="grid"/>.</param>
        public SimulationData(BaseBlock[] grid)
        {
            this.grid = grid;
        }
    }
}
