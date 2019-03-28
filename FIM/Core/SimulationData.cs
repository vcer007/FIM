using FIM.FluidData;
using FIM.Report;
using FIM.RockData;
using FIM.Well;

namespace FIM.Core
{

    /// <summary>
    /// This class contains all the general data of the simulation run "Run Specifications".
    /// </summary>
    public class SimulationData
    {
        // RunSpec

        /// <summary>
        /// The title of the simulation run.
        /// </summary>
        public string title;

        /// <summary>
        /// The solution procedure used for the model.
        /// </summary>
        /// <seealso cref="Global.SolutionProcedure"/>
        public Global.SolutionProcedure solutionProcedure;

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
        /// A boolean to indicate if vaporized is present.
        /// </summary>
        /// <seealso cref="phases"/>
        public bool vaporizedOilPresent;

        /// <summary>
        /// This variable determines if gas is allowed to re-dissolve in oil. This is equivalent to eclipse having no DRSDT keyword.
        /// </summary>
        public bool isGasResolutionAllowed;

        /// <summary>
        /// Model grid dimensions.
        /// </summary>
        public int x, y, z;

        // Grid

        // the grid "a list of the blocks"        
        /// <summary>
        /// An array of <see cref="BaseBlock"/> that compose the grid.
        /// </summary>
        public BaseBlock[] grid;

        // Props

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
        public PorosityCalculator porosityCalculator;

        // Solution
        // Summary        
        /// <summary>
        /// An instance of the <see cref="Output"/> class.
        /// </summary>
        public Output output;

        // Schedule

        // the grid "a list of the wells".     
        /// <summary>
        /// An array of <see cref="BaseWell"/> present in the model.
        /// </summary>
        public BaseWell[] wells;

        // Internal Simulator Configurations

        /// <summary>
        /// The original relaxation factor value.
        /// </summary>
        /// <remarks>
        /// The <see cref="relaxationFactor"/> is set to this value at the beginning of the time step.
        /// </remarks>
        /// <seealso cref="SimulationData.minimumRelaxation"/>
        public double originalRelaxationFactor;

        /// <summary>
        /// The minimum relaxation factor value.
        /// </summary>
        /// <remarks>
        /// <para>Default value is equal to 0.5.</para>
        /// <para>If the solution does not converge till this value, the time step is reduced.</para>
        /// </remarks>
        /// <seealso cref="SimulationData.timeStepSlashingFactor"/>
        /// <seealso cref="SimulationData.relaxationFactor"/>
        public double minimumRelaxation;

        /// <summary>
        /// The relaxation factor used when the solution oscillates.
        /// </summary>
        /// <remarks>
        /// When oscillations are detected, this value is used to "dampen" the calculated dx.
        /// </remarks>
        /// <seealso cref="originalRelaxationFactor"/>
        /// <seealso cref="minimumRelaxation"/>
        public double relaxationFactor;

        /// <summary>
        /// This value is used to decrement the relaxation factor if oscillation is detected.
        /// </summary>
        public double relaxationFactorDecrement;

        /// <summary>
        /// the maximum value of the ratio of material balance errors between two subseuent non-linear iterations.
        /// </summary>
        public double maximumMaterialBalanceErrorRatio;

        /// <summary>
        /// The original_time_step
        /// </summary>
        public double originalTimeStep;

        /// <summary>
        /// The time step used in calculations.
        /// </summary>
        public double timeStep;

        /// <summary>
        /// The factor used to reduce the <see cref="SimulationData.timeStep"/> when the solver stagnates or oscillates.
        /// </summary>
        public double timeStepSlashingFactor;

        /// <summary>
        /// The total time of the simulation.
        /// </summary>
        public double endingTime;

        /// <summary>
        /// The maximum number of non-linear iterations used.
        /// </summary>
        /// <remarks>
        /// The maximum number of iterations used in Newton-Raphson non-linear iterations befor the solver slashes the time step.
        /// </remarks>
        /// <seealso cref="timeStep"/>
        /// <seealso cref="timeStepSlashingFactor"/>
        public double maximumNonLinearIterations;

        /// <summary>
        /// The material balance error between fluid originally and currently in-place.
        /// </summary>
        /// <seealso cref="MBE_Tolerance"/>
        public double MBE_Oil, MBE_Gas, MBE_Water;

        /// <summary>
        /// The tolerance used to determine convergence.
        /// </summary>
        /// <remarks>
        /// currently tolerance is set as an absolute variance between original fluid in place and fluid in place.
        /// This behavior is used to ensure very high accuracy.
        /// However, a percentage may be used to yield sufficiently accurate solutions with a lower number of iterations required.
        /// </remarks>
        public double MBE_Tolerance;

        internal bool Gravity;

        // The class constructor.
        /// <summary>
        /// Initializes a new instance of the <see cref="SimulationData"/> class.
        /// </summary>
        /// <param name="grid">The <see cref="grid"/>.</param>
        public SimulationData(BaseBlock[] grid)
        {
            this.grid = grid;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="SimulationData"/> class.
        /// </summary>
        /// <remarks>
        /// This is an empty constructor used in the parser.
        /// </remarks>
        public SimulationData()
        {
        }
    }
}