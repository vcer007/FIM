using System;
using FIM.Core;

namespace FIM.MaterialBalance
{
    /// <summary>
    /// This class handles the calculation of material balance errors for the different phases.
    /// </summary>
    public static class MBE
    {

        /// <summary>
        /// calculates the Oil material balance error.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <returns>The value of MBE for the oil phase</returns>
        /// <seealso cref="Global.Phase"/>
        public static double CheckOil(SimulationData data)
        {
            double OOIP = 0, OIP = 0, q = 0;
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                OOIP += block.Vp[0] * block.So[0] / block.Bo[0];
                OIP += block.Vp[1] * block.So[1] / block.Bo[1];
            }

            OOIP += GetVaporizedOilInPlace(data, 0);
            OIP += GetVaporizedOilInPlace(data, 1);

            for (int i = 0; i < data.wells.Length; i++)
            {
                q += Global.a * (data.wells[i].q_oil[1] + data.wells[i].q_vap_oil[1]) * data.timeStep;
            }

            double difference = OOIP - (OIP + q);
            var pore_volume = GetTotalPoreVolume(1, data);

            return difference / /*(OIP + q) * percentage_factor*/ pore_volume;
        }

        private static double GetVaporizedOilInPlace(SimulationData data, int time_level)
        {
            double temp = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                temp += data.grid[i].Rvo[time_level] * data.grid[i].Vp[time_level] * data.grid[i].Sg[time_level] / data.grid[i].Bg[time_level];
            }

            return temp;
        }

        /// <summary>
        /// calculates the Gas material balance error.
        /// </summary>
        /// <remarks>
        /// This method calculate the material balance error for the whole gas phase.
        /// Gas soluble in oil is taken into consideration.
        /// </remarks>
        /// <param name="data">The data.</param>
        /// <returns>The value of MBE for the gas phase</returns>
        /// <seealso cref="Global.Phase"/>
        /// <seealso cref="SimulationData.solubleGasPresent"/>
        public static double CheckGas(SimulationData data)
        {
            double OGIP = 0, GIP = 0, q = 0;

            // original gas in place.
            OGIP = GetGIP(data, 0);

            // gas currently in place.
            GIP = GetGIP(data, 1);

            for (int i = 0; i < data.wells.Length; i++)
            {
                q += Global.a * (data.wells[i].q_free_gas[1] + data.wells[i].q_solution_gas[1]);
            }

            double difference = OGIP - (GIP + q * data.timeStep);
            var pore_volume = GetTotalPoreVolume(1, data);

            return difference / /*(GIP + q * data.timeStep) * percentage_factor*/pore_volume;
        }

        /// <summary>
        /// calculates the Water material balance error.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <returns>The value of MBE for the water phase</returns>
        /// <seealso cref="Global.Phase"/>
        public static double CheckWater(SimulationData data)
        {
            double OWIP = 0, WIP = 0, q = 0;
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                OWIP += block.Vp[0] * block.Sw[0] / block.Bw[0];
                WIP += block.Vp[1] * block.Sw[1] / block.Bw[1];
            }

            for (int i = 0; i < data.wells.Length; i++)
            {
                q += Global.a * data.wells[i].q_water[1] * data.timeStep;
            }

            double difference = OWIP - (WIP + q);
            var pore_volume = GetTotalPoreVolume(1, data);

            return difference  / /*(WIP + q) * percentage_factor*/pore_volume;
        }

        /// <summary>
        /// calculates the gas in place at a given time level in units of cubic foot.
        /// </summary>
        /// <param name="data"><see cref="SimulationData"/></param>
        /// <param name="timeLevel">The time level needed; n0 or n1 "0 or 1"</param>
        /// <returns>The value of GIP</returns>
        public static double GetGIP(SimulationData data, int timeLevel)
        {
            double temp = FreeGasIP(data, timeLevel) + SolubleGasIP(data, timeLevel);

            return temp;
        }

        /// <summary>
        /// calculates the free gas in place at a given time level in units of cubic foot.
        /// </summary>
        /// <param name="data"><see cref="SimulationData"/></param>
        /// <param name="timeLevel">The time level needed; n0 or n1 "0 or 1"</param>
        /// <returns>The value of FreeGasIP</returns>
        public static double FreeGasIP(SimulationData data, int timeLevel)
        {
            double temp = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                temp += data.grid[i].Vp[timeLevel] * data.grid[i].Sg[timeLevel] / data.grid[i].Bg[timeLevel];
            }

            return temp;
        }

        /// <summary>
        /// calculates the soluble gas in place at a given time level in units of cubic foot.
        /// </summary>
        /// <param name="data"><see cref="SimulationData"/></param>
        /// <param name="timeLevel">The time level needed; n0 or n1 "0 or 1"</param>
        /// <returns>The value of SolubleGasIP</returns>
        public static double SolubleGasIP(SimulationData data, int timeLevel)
        {
            double temp = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                temp += data.grid[i].Rso[timeLevel] * data.grid[i].Vp[timeLevel] * data.grid[i].So[timeLevel] / data.grid[i].Bo[timeLevel];
            }

            return temp;
        }

        private static double GetTotalPoreVolume(int time_level, SimulationData data)
        {
            double poreVolume = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                poreVolume += data.grid[i].Vp[time_level];
            }

            return poreVolume;
        }
    }
}
