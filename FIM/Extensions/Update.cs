using FIM.Core;

namespace FIM.Extensions
{
    /// <summary>
    /// An extension class that contains extension methods to <see cref="BaseBlock"/> for updating the properties.
    /// </summary>
    /// <remarks>
    /// <para>The extension methods in this class are used for updating the <see cref="BaseBlock"/> member variables array in accordance to <see cref="Global.STEPS_MEMORY"/>.</para>
    /// <para>The update methods have been implemented as extensions rather than internal member methods in the <see cref="BaseBlock"/> class to reduce clutter and to allow for flexibility in adding new methods in a similar way.</para>
    /// </remarks>
    public static class Update
    {
        /// <summary>
        /// Updates the properties of the <see cref="BaseBlock"/>.
        /// </summary>
        /// <remarks>
        /// <para>Different time levels properties are updated recursively down-up.</para>
        /// <para>If time level is equal to zero, then n0, n1 and n2, if available, time levels will be updated.</para>
        /// <para>If time level is equal to 1, then only n1 and n2, if available, time levels will be updated.</para>
        /// </remarks>
        /// <param name="block">
        /// This is not input to the method. This is only used to efine the extension method
        /// Extension methods are used exactly in the same way as classes member methods.
        /// </param>
        /// <param name="data">The <seealso cref="SimulationData"/> object thaat contains all the input data.</param>
        /// <param name="P">The new pressure value.</param>
        /// <param name="Sw">The new water saturation value.</param>
        /// <param name="Sg">The new gas saturation value.</param>
        /// <param name="So">The new oil saturation value.</param>
        /// <param name="time_level"><para>Default value is 0.</para>
        /// <para>The time_level at which the properties should be updated.</para>
        /// </param>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        /// <seealso cref="Reset(BaseBlock, SimulationData)"/>
        public static void UpdateProperties(this BaseBlock block, SimulationData data, double P, double Sw, double Sg, double So, int time_level = 0)
        {
            // pressure
            block.P_previousStep = block.P[0];

            block.P[time_level] = P;

            // data.porosity_calculator
            block.porosity[time_level] = data.porosityCalculator.getPorosity(P);

            var new_rs = data.pvt.GetRs(Global.Phase.Oil, P);
            block.Rso[time_level] = new_rs;
            //if (time_level == 2)
            //{
            //    block.Rso[2] = block.Rso[1];
            //}
            //else
            //{
            //    if (block.Sg[time_level] == 0 || new_rs < block.Rso[0])
            //    {
            //        block.Rso[time_level] = new_rs;
            //    }
            //}

            block.So[time_level] = So;
            block.Sw[time_level] = Sw;
            block.Sg[time_level] = Sg;

            // update capillary first
            block.Pcgo[time_level] = data.pvt.GetGasCapillaryPressure(Sg);
            //block.Pcow[time_level] = data.pvt.GetWaterCapillaryPressure(Sw);

            // fluid
            block.Bo[time_level] = data.pvt.GetFVF(Global.Phase.Oil, P);
            block.Bw[time_level] = data.pvt.GetFVF(Global.Phase.Water, block.GetPw(time_level, time_level));
            block.Bg[time_level] = data.pvt.GetFVF(Global.Phase.Gas, block.GetPg(time_level, time_level));

            block.viscosityOil[time_level] = data.pvt.GetViscosity(Global.Phase.Oil, P);
            block.viscosityWater[time_level] = data.pvt.GetViscosity(Global.Phase.Water, /*block.GetPw(time_level)*/P);
            block.viscosityGas[time_level] = data.pvt.GetViscosity(Global.Phase.Gas, block.GetPg(time_level, time_level));

            //// Kro is only dependent on Sg.
            //block.Kro[time_level] = data.scal.GetKr(Global.Phase.Oil, Sg);
            //// Krw is dependent on sw
            //block.Krw[time_level] = data.scal.GetKr(Global.Phase.Water, Sw);
            //block.Krg[time_level] = data.scal.GetKr(Global.Phase.Gas, Sg);

            // Kro is only dependent on Sg.
            block.Kro[time_level] = data.scal.GetKro(Sg, Sw, 0.12);
            // Krw is dependent on sw
            block.Krw[time_level] = data.scal.GetKrw(Sw);
            block.Krg[time_level] = data.scal.GetKrg(Sg);

            // volumetric
            block.Vp[time_level] = block.bulkVolume * block.porosity[time_level];

            // update wells that are located within the block, if any.
            if (block.type == Global.BlockType.WellBlock)
            {
                for (int i = 0; i < data.wells.Length; i++)
                {
                    if (data.wells[i].index == block.index)
                    {
                        data.wells[i].Update(time_level, block, time_level == 2 ? true : false);
                    }
                }
            }

            if (time_level == 0)
            {
                // block way we update n1 time level properties also automatically
                // after updating n0 time level properties.
                // note that the initial guess here for n1 is assumed to be the same values as the ones at n0.
                block.UpdateProperties(data, P, Sw, Sg, So, 1);
            }
            else if (time_level == 1)
            {
                // check if it is a fully implicit simulation.
                // So we can update three time levels storage "where the third one, n2, is used solely for perturbation".
                if (data.solutionProcedure == Global.SolutionProcedure.FullyImplicit)
                {
                    block.UpdateProperties(data, P + Global.EPSILON_P, Sw + Global.EPSILON_S, Sg + Global.EPSILON_S, So, 2);
                }
            }
        }

        /// <summary>
        /// Resets n1, and n2 time level if available, properties to n0.
        /// </summary>
        /// <param name="block">
        /// This is not input to the method. This is only used to efine the extension method
        /// Extension methods are used exactly in the same way as classes member methods.
        /// </param>
        /// <param name="data">The data.</param>
        /// <seealso cref="SimulationData"/>
        /// <seealso cref="SimulationData.timeStepSlashingFactor"/>
        public static void Reset(this BaseBlock block, SimulationData data)
        {
            // to reset the properties of the block, simply call updateProperties_n0_n1 with input parameters from
            // the n0 time level.
            block.UpdateProperties(data, block.P[0], block.Sw[0], block.Sg[0], block.So[0], 0);
        }
    }
}
