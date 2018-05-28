using System;

namespace FIM.Core
{
    /// <summary>
    /// This class contains all constants and enumerations declarations.
    /// </summary>
    public static class Global
    {
        // Constants.

        /// <summary>
        /// Conversion factor used in the darcy transmissibility equation to get flow rates in unit bbl. per day.
        /// </summary>
        public const double Bc = 0.001127;

        /// <summary>
        /// Conversion factor from bbl. to cubic foot.
        /// </summary>
        public const double a = 5.614583;

        /// <summary>
        /// The epsilon used for caluclating derivatives numerically.
        /// </summary>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public static double EPSILON_P;
        public static double EPSILON_S;
        public static double MINIMUM;

        /// <summary>
        /// The size of the properties <see cref="Array"/>.
        /// </summary>
        /// <remarks>
        /// <para>Instead of a single value for each property "e.g; oil saturation", an array is used 
        /// to store different values for each time step.</para>
        /// <para>Typically, the time steps are n0, n+1 and n+1.</para>
        /// <para>The n+2 time level is used for perturbation.</para>
        /// <para>Independent variables "e.g; pressure, saturation, ..." will have their n+2 values set to n+1 + <see cref="Global.EPSILON_P"/></para>
        /// <para>Dependent variables "e.g; Bo, Bg, Kr, ..." will have their n+2 values set to
        /// corresponding values to their corresponding independent variable at the n+2 time level.</para>
        /// <para>This way adds more flexibility in using values at different time levels whenever needed.</para>
        /// <para>So, the typical array size will be 3.</para>
        /// </remarks>
        public const int STEPS_MEMORY = 3;

        // Enumerations.

        /// <summary>
        /// The different phases available for a black-oil system.
        /// </summary>
        public enum Phase {Water, Oil, Gas/*, DissolvedGas*/}

        /// <summary>
        /// The different types of a well according to production, or injection, status.
        /// </summary>
        public enum WellType {Injection, Production, ShutIn}

        /// <summary>
        /// The different types of well control schemes.
        /// </summary>
        public enum WellControl { OilRate, GasRate, WaterRate, BHP}

        /// <summary>
        /// The method used in calculating well flow rates.
        /// </summary>
        /// <remarks>
        /// Implicit means that the time level of the block pressure used for calculating the drawdown is n1.
        /// Explicit means using n0 time level of the block pressure for calculating the well drawdown.
        /// </remarks>
        public enum WellRateCalculation { /** Explicit Calculation*/Explicit, /** Implicit calculation*/Implicit}

        /// <summary>
        /// Indicates whether the block contains a well or not.
        /// </summary>
        public enum BlockType { /** A well Block*/WellBlock, /** An ordinary block*/NormalBlock, /** An inactive Block*/InactiveBlock }

        /// <summary>
        /// Sections in eclipse input data file.
        /// </summary>
        public enum DataFileSection { /** The Run Specifications*/RunSpec, /** The grid*/Grid, /** The fluids properties*/Props, /** Initial equilibrium conditions*/Solution, /** Repoerts*/Summary, /** Well data*/Schedule }

        /// <summary>
        /// The solution methods available.
        /// </summary>
        public enum SolutionProcedure { /** A fully implicit solver*/FullyImplicit, /** An IMPES solver*/IMPES }

        /// <summary>
        /// Units used for the data supplied in input data files.
        /// </summary>
        public enum UnitConvention { /** Field unit system*/Field, /** Metric unit system*/Metric, /** Lab unit system*/Lab }

        /// <summary>
        /// Methods for interpolation.
        /// </summary>
        public enum Interpolation {/** Table look up and interpolation*/TableLookUp, /** Using a predefined polynomial*/PolynomialFunction, /** Using a predefined power function*/PowerFunction }

        /// <summary>
        /// The independent variables accounted for in numerical perturbation.
        /// </summary>
        public enum Variable {/** Derivative with respect to pressure*/Pressure, /** Derivative with respect to gas saturation*/SaturationGas, /** Derivative with respect to water saturation*/SaturationWater }

        /// <summary>
        /// The padding used to format the output text.
        /// </summary>
        public static int padding;

        /// <summary>
        /// The number of the decimal places of the output values.
        /// </summary>
        public static string decimalPlaces;

        public static double gamma_c = 0.21584E-3;

        public static double alpha = 32.174;
    }
}
