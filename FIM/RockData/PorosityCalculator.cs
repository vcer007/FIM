namespace FIM.RockData
{
    /// <summary>
    /// This class is used to calculate porosity_calculator porosities at different pressures.
    /// </summary>
    /// <remarks>
    /// <para>The <see cref="porosity_ref"/> member variable in this class is used to store porosity_calculator at reference pressure <see cref="pressure_ref"/>.</para>
    /// <para><see cref="Cf"/> is the porosity_calculator compressibility used in the correlation to get porosities at different pressures.</para>
    /// </remarks>
    /// <seealso cref="FluidData.PVT"/>
    /// <seealso cref="FluidData.SCAL"/>
    public class PorosityCalculator
    {
        //Rock compressibility factor
        double Cf;

        //Value of the porosity at the reference pressure.
        double porosity_ref;

        //The reference pressure at which the porosity was measured.
        double pressure_ref;

        /// <summary>
        /// Initializes a new instance of the <see cref="PorosityCalculator"/> class.
        /// </summary>
        /// <param name="Cf">The porosity_calculator compressibility.</param>
        /// <param name="porosity_ref">The porosity_calculator value at reference pressure.</param>
        /// <param name="pressure_ref">The reference pressure at which the reference porosity_calculator was measured.</param>
        public PorosityCalculator(double Cf, double porosity_ref, double pressure_ref)
        {
            this.Cf = Cf;
            this.porosity_ref = porosity_ref;
            this.pressure_ref = pressure_ref;
        }

        /// <summary>
        /// Gets the porosity_calculator at a certain value of pressure.
        /// </summary>
        /// <param name="pressure">The new pressure.</param>
        /// <returns>The new value of porosity_calculator at the newly specified pressure.</returns>
        public double getPorosity(double pressure)
        {
            return porosity_ref * (1 + Cf * (pressure - pressure_ref));
        }

    }
}
