using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.MaterialBalance
{
    class MBE
    {
        /// <summary>
        /// calculates the Oil material balance error.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <returns>The value of MBE for the oil phase</returns>
        /// <seealso cref="Global.Phase"/>
        public static double checkOil(SimulationData data)
        {
            double tolerance = 1e-7;

            double OOIP = 0, OIP = 0, q = 0;
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                OOIP += block.Vp[0] * block.So[0] / block.Bo[0];
                OIP += block.Vp[1] * block.So[1] / block.Bo[1];

                //if (block.type == Global.BlockType.Well_Block)
                //{
                //    q += Global.a * block.q_oil[1] * data.time_step;
                //}
            }

            for (int i = 0; i < data.wells.Length; i++)
            {
                q += Global.a * data.wells[i].q_oil[1] * data.time_step;
            }

            return OOIP - (OIP + q);
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
        public static double checkGas(SimulationData data)
        {
            double tolerance = 1e-7;

            double OGIP = 0, GIP = 0, q = 0;
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                OGIP += block.Vp[0] * block.Sg[0] / block.Bg[0];
                GIP += block.Vp[1] * block.Sg[1] / block.Bg[1];

                //if (block.type == Global.BlockType.Well_Block)
                //{
                //    if (block.well_type == Global.WellType.Production)
                //    {
                //        q += Global.a * block.q_gas[1];
                //    }
                //    else
                //    {
                //        q += Global.a * block.specified_flow_rate;
                //    }
                //}

                if (data.solubleGasPresent)
                {
                    OGIP += block.Rso[0] * block.Vp[0] * block.So[0] / block.Bo[0];
                    GIP += block.Rso[1] * block.Vp[1] * block.So[1] / block.Bo[1];

                    //q += Global.a * block.Rso[1] * block.q_oil[1];
                    //q += Global.a * block.Rso[1] * block.q_oil[1];
                }
            }

            for (int i = 0; i < data.wells.Length; i++)
            {
                q += Global.a * (data.wells[i].q_free_gas[1] + data.wells[i].q_solution_gas[1]);
            }

            //q = 0;
            return OGIP - (GIP + q * data.time_step);
        }

        /// <summary>
        /// calculates the Water material balance error.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <returns>The value of MBE for the water phase</returns>
        /// <seealso cref="Global.Phase"/>
        internal static double checkWater(SimulationData data)
        {
            double tolerance = 1e-7;

            double OWIP = 0, WIP = 0, q = 0;
            BaseBlock block;

            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                OWIP += block.Vp[0] * block.Sw[0] / block.Bw[0];
                WIP += block.Vp[1] * block.Sw[1] / block.Bw[1];

                //if (block.type == Global.BlockType.Well_Block)
                //{
                //    q += Global.a * block.q_water[1] * data.time_step;
                //}
            }

            for (int i = 0; i < data.wells.Length; i++)
            {
                q += Global.a * data.wells[i].q_water[1] * data.time_step;
            }

            return OWIP - (WIP + q);
        }

        public static double GIP(SimulationData data, int time_level)
        {
            BaseBlock block;
            double temp = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                temp += block.Vp[time_level] * block.Sg[time_level] / block.Bg[time_level];

                if (data.solubleGasPresent)
                {
                    temp += block.Rso[time_level] * block.Vp[time_level] * block.So[time_level] / block.Bo[time_level];
                }
            }

            return temp;
        }

        public static double FreeGasIP(SimulationData data)
        {
            double temp = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                temp += data.grid[i].Vp[1] * data.grid[i].Sg[1] / data.grid[i].Bg[1];
            }

            return temp;
        }

        public static double SolubleGasIP(SimulationData data)
        {
            double temp = 0;

            for (int i = 0; i < data.grid.Length; i++)
            {
                temp += data.grid[i].Rso[1] * data.grid[i].Vp[1] * data.grid[i].So[1] / data.grid[i].Bo[1];
            }

            return temp;
        }
    }
}
