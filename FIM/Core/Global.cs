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
    public class Global
    {
        // Constants.

        public const double Bc = 0.001127;
        public const double a = 5.614583;
        //public const double a = 5.615;
        public const double epsilon = 1E-5;
        public const double initial_Sg = 0;
        public const double epsilon_p = 1E-5;
        public const double epsilon_s = 1E-5;
        public const double PI = Math.PI;

        // Enumerations.

        public enum Phase {Water, Oil, Gas, SolubleGas}
        public enum WellType {Injection, Production, ShutIn}
        public enum BlockType { Well_Block, Normal_Block, Inactive_Block }
        public enum DataFileSection { RunSpec }
        public enum SolutionProcedure { IMPES, Fully_Implicit }
        public enum UnitConvention { Field, Metric, Lab }
        public enum LinearInterpolation {TableLookUp, PolynomialFunction, PowerFunction}
        public enum Variable {Pressure, Saturation_Oil, Saturation_Gas, Saturation_Water}

        // Delegates "functions with a defined form and (return and input) types but with no specific form"

        // The general form of the interpolation equation used to calculate PVT and rock properties from pressure "or saturation".
        // "input" is pressure "for FVF, viscosity and Vp" and saturation "for Kr".
        public delegate double getProperty(double input);
    }
}
