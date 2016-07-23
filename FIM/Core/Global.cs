﻿using System;

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
        public const double EPSILON = 1E-5;

        /// <summary>
        /// The size of the properties <see cref="Array"/>.
        /// </summary>
        /// <remarks>
        /// <para>Instead of a single value for each property "e.g; oil saturation", an array is used 
        /// to store different values for each time step.</para>
        /// <para>Typically, the time steps are n0, n+1 and n+1.</para>
        /// <para>The n+2 time level is used for perturbation.</para>
        /// <para>Independent variables "e.g; pressure, saturation, ..." will have their n+2 values set to n+1 + <see cref="Global.EPSILON"/></para>
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
        public enum Phase {Water, Oil, Gas, SolubleGas}

        /// <summary>
        /// The different types of a well according to production, or injection, status.
        /// </summary>
        public enum WellType {Injection, Production, ShutIn}

        /// <summary>
        /// The different types of well control schemes.
        /// </summary>
        public enum WellControl { OilRate, GasRate, BHP}

        /// <summary>
        /// Indicates whether the block contains a well or not.
        /// </summary>
        public enum BlockType { WellBlock, NormalBlock, InactiveBlock }

        /// <summary>
        /// Sections in eclipse input data file.
        /// </summary>
        public enum DataFileSection { RunSpec, Grid, Props, Solution, Summary, Schedule }

        /// <summary>
        /// The solution methods available.
        /// </summary>
        public enum SolutionProcedure { FullyImplicit, IMPES }

        /// <summary>
        /// Units used for the data supplied in input data files.
        /// </summary>
        public enum UnitConvention { Field, Metric, Lab }

        /// <summary>
        /// Methods for interpolation.
        /// </summary>
        public enum Interpolation {TableLookUp, PolynomialFunction, PowerFunction}

        /// <summary>
        /// The independent variables accounted for in numerical perturbation.
        /// </summary>
        public enum Variable {Pressure, SaturationOil, SaturationGas, SaturationWater}
    }
}
