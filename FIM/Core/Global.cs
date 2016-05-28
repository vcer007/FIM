using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Core
{
    /// <summary>
    /// This class contains all constants and enumerations declarations.
    /// </summary>
    class Global
    {
        // Constants.

        public const double Bc = 0.001127;
        public const double a = 5.614583;
        public const double epsilon = 0.0000001;

        // Enumerations.

        public enum Phase {Water, Oil, Gas}
        public enum WellType {Injection, Production, ShutIn}
        public enum BlockType { Well_Block, Normal_Block, Inactive_Block }
        public enum DataFileSection { RunSpec }
        public enum SolutionProcedure { IMPES, Fully_Implicit }
        public enum UnitConvention { Field, Metric, Lab }
        public enum LinearInterpolation {TableLookUp, PolynomialFunction, PowerFunction}
        public enum Variable {Pressure, Saturation}

        // Delegates "functions with a defined form and (return and input) types but with no specific form"

        // The general form of the interpolation equation used to calculate PVT and rock properties from pressure "or saturation".
        // "input" is pressure "for FVF, viscosity and Vp" and saturation "for Kr".
        public delegate double getProperty(double input);
    }
}
