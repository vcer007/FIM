using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Extensions
{
    public static class SingleBlockExtensions
    {
        private static double COP, COG, CGP, CGG, q_oil, q_free_gas, q_soluble_gas, temp, bhp;

        public static double getCOP(this BaseBlock block, SimulationData data)
        {
            COP = 1 / (Global.a * data.time_step) * (block.phi_dash() / block.Bo[0] + block.Vp[1] * block.FVF_dash(Global.Phase.Oil)) * (1 - block.Sg[0]);

            //q_oil = Well.WellData.calculateFlow_Rate(block.P[2], block.BHP[1], block.Kro[1], block.viscosity_oil[2], block.WI, block.Bo[2]);
            //Console.WriteLine(q_oil);

            //temp = (q_oil - block.q_oil[1]) / Global.epsilon_p;

            return - COP;
        }
        public static double getCOG(this BaseBlock block, SimulationData data)
        {
            COG = -1 / (Global.a * data.time_step) * (block.Vp[1] / block.Bo[1]);

            //q_oil = Well.WellData.calculateFlow_Rate(block.P[1], block.BHP[1], block.Kro[2], block.viscosity_oil[1], block.WI, block.Bo[1]);

            //temp = (q_oil - block.q_oil[1]) / Global.epsilon_s;

            return - COG;
        }
        public static double getCGP(this BaseBlock block, SimulationData data)
        {
            CGP = 1 / (Global.a * data.time_step) * (((block.phi_dash() / block.Bo[0] + block.Vp[1] * block.FVF_dash(Global.Phase.Oil)) * block.Rso[0] + (block.Vp[1] / block.Bo[1]) * block.Rso_dash()) * (1 - block.Sg[0]) + (block.phi_dash() / block.Bg[0] + block.Vp[1] * block.FVF_dash(Global.Phase.Gas)) * block.Sg[0]);

            //q_soluble_gas = block.Rso[2] * Well.WellData.calculateFlow_Rate(block.P[2], block.BHP[1], block.Kro[1], block.viscosity_oil[2], block.WI, block.Bo[2]);
            //temp = (q_soluble_gas - block.Rso[1] * block.q_oil[1]) / Global.epsilon_p;

            temp = (block.Rso[2] - block.Rso[1]) / (block.P[2] - block.P[1]) * block.q_oil[1];

            bhp = Well.WellData.calculatePwf(block, block.P[2], block.Kro[1], block.viscosity_oil[2], block.Bo[2]);
            temp += (Well.WellData.calculateFlow_Rate(block.P[2], bhp, block.Krg[1], block.viscosity_gas[2], block.WI, block.Bg[2]) - block.q_gas[1]) / (block.P[2] - block.P[1]);

            return - CGP - temp;
        }
        public static double getCGG(this BaseBlock block, SimulationData data)
        {
            CGG = 1 / (Global.a * data.time_step) * (block.Vp[1] / block.Bg[1] - (block.Vp[1] / block.Bo[1]) * block.Rso[1]);

            //q_soluble_gas = block.Rso[1] * Well.WellData.calculateFlow_Rate(block.P[1], block.BHP[1], block.Kro[2], block.viscosity_gas[1], block.WI, block.Bg[1]);
            //temp = (q_soluble_gas - block.Rso[1] * block.q_oil[1]) / Global.epsilon_s;
            //temp = 0;
            bhp = Well.WellData.calculatePwf(block, block.P[1], block.Kro[2], block.viscosity_oil[1], block.Bo[1]);
            temp = (Well.WellData.calculateFlow_Rate(block.P[1], bhp, block.Krg[2], block.viscosity_gas[1], block.WI, block.Bg[1]) - block.q_gas[1]) / (block.Sg[2] - block.Sg[1]);

            return - CGG - temp;
        }

        //public static double getQ_P(this BaseBlock block, SimulationData data)
        //{

        //}

        public static double getminus_R(this BaseBlock block, SimulationData data, Global.Phase phase)
        {
            double R = 0;

            if (phase == Global.Phase.Oil)
            {
                R = -1 / (Global.a * data.time_step) * (block.Vp[1] * (1 - block.Sg[1]) / block.Bo[1] - block.Vp[0] * (1 - block.Sg[0]) / block.Bo[0]) - block.q_oil[1];
            }
            else if (phase == Global.Phase.Gas)
            {
                R = -1 / (Global.a * data.time_step) * (block.Vp[1] * block.Sg[1] / block.Bg[1] - block.Vp[0] * block.Sg[0] / block.Bg[0]) - block.q_gas[1];

                if (data.solubleGasPresent)
                {
                    R += -1 / (Global.a * data.time_step) * (block.Rso[1] * block.Vp[1] * (1 - block.Sg[1]) / block.Bo[1] - block.Rso[0] * block.Vp[0] * (1 - block.Sg[0]) / block.Bo[0]) - block.Rso[1] * block.q_oil[1];
                }
            }

            return -1 * R;
        }

        private static double phi_dash(this BaseBlock block)
        {
            double P_difference = block.P[1] - block.P[0];
            return (block.Vp[1] - block.Vp[0]) / (P_difference);
        }

        private static double Rso_dash(this BaseBlock block)
        {
            double P_difference = block.P[1] - block.P[0];
            return (block.Rso[1] - block.Rso[0]) / (P_difference);
        }

        private static double FVF_dash(this BaseBlock block, Global.Phase phase)
        {
            double P_difference = block.P[1] - block.P[0];

            double FVF_1 = 1, FVF_0 = 1;
            if (phase == Global.Phase.Oil)
            {
                FVF_1 = block.Bo[1];
                FVF_0 = block.Bo[0];
            }
            else if (phase == Global.Phase.Water)
            {
                FVF_1 = block.Bw[1];
                FVF_0 = block.Bw[0];
            }
            else
            {
                FVF_1 = block.Bg[1];
                FVF_0 = block.Bg[0];
            }
            return (1 / FVF_1 - 1 / FVF_0) / (P_difference);
        }

        public static double[] calculate_minus_R(SimulationData data)
        {
            double[] minus_R = new double[data.grid.Length * data.phases.Length];

            BaseBlock block;
            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                minus_R[counter] = block.getminus_R(data, Global.Phase.Oil);
                minus_R[counter + 1] = block.getminus_R(data, Global.Phase.Gas);

                counter += data.phases.Length;
            }

            return minus_R;
        }

        public static double[][] calculateJacobians(SimulationData data, double[] minus_R)
        {
            int size = data.x * data.y * data.z * data.phases.Length;
            double[][] jacobians = new double[size][];

            BaseBlock block;

            int counter = 0;
            for (int i = 0; i < data.grid.Length; i++)
            {
                block = data.grid[i];

                jacobians[counter] = new double[size];
                jacobians[counter + 1] = new double[size];

                #region Oil
                // with respect to P
                jacobians[counter][data.phases.Length * block.index] = block.getCOP(data);
                // with respect to Sg
                jacobians[counter][data.phases.Length * block.index + 1] = block.getCOG(data);
                // with respect to Sw
                #endregion
                #region Gas
                // with respect to P
                jacobians[counter + 1][data.phases.Length * block.index] = block.getCGP(data);
                // with respect to Sg
                jacobians[counter + 1][data.phases.Length * block.index + 1] = block.getCGG(data);

                #endregion

                counter += data.phases.Length;
            }

            return jacobians;

        }

    }
}
