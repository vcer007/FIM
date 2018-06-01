using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Well
{
    /// <summary>
    /// This class contains the definition of a well's properties and associated methods.
    /// </summary>
    public class BaseWell
    {
        /// <summary>
        /// The name of the well.
        /// </summary>
        public string name;

        /// <summary>
        /// The index of the block in which the well is located.
        /// </summary>
        public int index;

        /// <summary>
        /// The radius of the well.
        /// </summary>
        public double radius;

        /// <summary>
        /// The skin factor.
        /// </summary>
        public double skin;

        /// <summary>
        /// The type of the well.
        /// </summary>
        /// <seealso cref="Global.WellType"/>
        public Global.WellType type;

        /// <summary>
        /// The well control method.
        /// </summary>
        /// <seealso cref="Global.WellControl"/>
        public Global.WellControl control;

        /// <summary>
        /// The method used for calculating well rates.
        /// </summary>
        /// <seealso cref="Global.WellRateCalculation"/>
        public Global.WellRateCalculation method;

        /// <summary>
        /// The specified minimum BHP
        /// </summary>
        public double specifiedMinimumBHP;

        /// <summary>
        /// The specified flow rate
        /// </summary>
        public double specifiedFlowRate;

        /// <summary>
        /// The equivalent well radius used for calculating flow rates in peaceman's well model.
        /// </summary>
        public double R_equivalent;

        /// <summary>
        /// The well productivity index.
        /// </summary>
        public double WI;

        /// <summary>
        /// The bottom hole pressure of the well.
        /// </summary>
        public double[] BHP;

        /// <summary>
        /// The well flow rates.
        /// </summary>
        public double[] q_oil, q_vap_oil, q_free_gas, q_solution_gas, q_water;

        // derivatives

        /// <summary>
        /// The oil flow rates derivatives.
        /// </summary>
        public double dq_oil_dP, dq_oil_dSg, dq_oil_dSw;

        public double dq_vap_oil_dP, dq_vap_oil_dSg, dq_vap_oil_dSw;

        /// <summary>
        /// The free gas flow rates derivatives.
        /// </summary>
        public double dq_free_gas_dP, dq_free_gas_dSg, dq_free_gas_dSw;

        /// <summary>
        /// The solution gas flow rates derivatives.
        /// </summary>
        public double dq_solution_gas_dP, dq_solution_gas_dSg, dq_solution_gas_dSw;

        /// <summary>
        /// The water flow rates derivatives.
        /// </summary>
        public double dq_water_dP, dq_water_dSg, dq_water_dSw;

        public SimulationData data { get; set; }

        /// <summary>
        /// Initializes a new instance of the <see cref="BaseWell"/> class.
        /// </summary>
        /// <param name="data">The <see cref="SimulationData"/>.</param>
        /// <param name="index">The index of the block the well is located into.</param>
        /// <param name="type">The well type<see cref="Global.WellType"/>.</param>
        /// <param name="control">The well control method <see cref="Global.WellControl"/>.</param>
        /// <param name="radius">The well radius.</param>
        /// <param name="skin">The skin factor.</param>
        /// <param name="specifiedMinimumBHP">The specified_minimum BHP.</param>
        /// <param name="specifiedFlowRate">The specified flow rate.</param>
        /// <param name="method">The well rate method used <see cref="Global.WellRateCalculation"/>.</param>
        public BaseWell(SimulationData data, int index, Global.WellType type, Global.WellControl control, double radius, double skin, double specifiedMinimumBHP = 0, double specifiedFlowRate = 0, Global.WellRateCalculation method = Global.WellRateCalculation.Explicit)
        {
            this.index = index;
            this.type = type;
            this.control = control;
            this.method = method;
            data.grid[index].type = Global.BlockType.WellBlock;
            this.radius = radius;
            this.skin = skin;
            this.R_equivalent = GetR_Equivalent(data.grid[index], data.grid);
            this.WI = GetWI(data.grid[index]);
            this.specifiedMinimumBHP = specifiedMinimumBHP;

            this.specifiedFlowRate = specifiedFlowRate;
            if (type == Global.WellType.Injection && control == Global.WellControl.GasRate)
            {
                this.specifiedFlowRate = -1 * specifiedFlowRate * 1000 / Global.a; // the input is in units of MSCF, but the internal units are unified to be bbl/day
            }
            else if (type == Global.WellType.Injection && control == Global.WellControl.WaterRate)
            {
                this.specifiedFlowRate = -1 * specifiedFlowRate;
            }

            this.BHP = new double[Global.STEPS_MEMORY];
            this.q_oil = new double[Global.STEPS_MEMORY];
            this.q_vap_oil = new double[Global.STEPS_MEMORY];
            this.q_free_gas = new double[Global.STEPS_MEMORY];
            this.q_solution_gas = new double[Global.STEPS_MEMORY];
            this.q_water = new double[Global.STEPS_MEMORY];

            this.data = data;
        }

        /// <summary>
        /// Updates the well properties at the specified time level.
        /// </summary>
        /// <param name="timeLevel">The time_level.</param>
        /// <param name="block">The block.</param>
        /// <param name="derivatives">if set to <c>true</c> [derivatives].</param>
        public void Update(int timeLevel, BaseBlock block, bool derivatives = true)
        {
            UpdateFlowRates(timeLevel, block);
            if (CheckWellControlChanged(timeLevel))
            {
                UpdateFlowRates(timeLevel, block);
            }

            //only for fully implicit simulation do we need to calculate well rates derivatives.
            if (derivatives)
            {
                UpdateDerivatives(block);
            }
        }

        // updates well flow rates at a certain time level.
        private void UpdateFlowRates(int timeLevel, BaseBlock block)
        {
            if (control == Global.WellControl.OilRate)
            {
                if (type == Global.WellType.Production)
                {
                    q_oil[timeLevel] = specifiedFlowRate;
                    BHP[timeLevel] = CalculatePwf(block.P[timeLevel], block.Kro[timeLevel], block.viscosityOil[timeLevel], block.Bo[timeLevel]);
                    q_free_gas[timeLevel] = CalculateFlowRate(block.GetPg(timeLevel, timeLevel), BHP[timeLevel], block.Krg[timeLevel], block.viscosityGas[timeLevel], WI, block.Bg[timeLevel]);
                    q_vap_oil[timeLevel] = block.Rvo[timeLevel] * q_free_gas[timeLevel];
                    q_solution_gas[timeLevel] = block.Rso[timeLevel] * q_oil[timeLevel];
                    q_water[timeLevel] = CalculateFlowRate(block.GetPw(timeLevel, timeLevel), BHP[timeLevel], block.Krw[timeLevel], block.viscosityWater[timeLevel], WI, block.Bw[timeLevel]);
                }
                else if (type == Global.WellType.Injection)
                {
                    //To-Do : 
                }
                else
                {
                    // shut-in well
                    q_oil[timeLevel] = 0;
                    q_vap_oil[timeLevel] = 0;
                    q_free_gas[timeLevel] = 0;
                    q_solution_gas[timeLevel] = 0;
                    q_water[timeLevel] = 0;
                }

            }
            else if (control == Global.WellControl.GasRate)
            {
                if (type == Global.WellType.Production)
                {
                    //To-Do : 
                }
                else if (type == Global.WellType.Injection)
                {
                    q_free_gas[timeLevel] = specifiedFlowRate;
                }
                else
                {
                    // shut-in well
                    //q_oil[time_level] = 0;
                    q_free_gas[timeLevel] = 0;
                    q_vap_oil[timeLevel] = 0;
                    //q_solution_gas[time_level] = 0;
                    //q_water[time_level] = 0;
                }
            }
            else if (control == Global.WellControl.WaterRate)
            {
                if (type == Global.WellType.Production)
                {
                    //To-Do : 
                }
                else if (type == Global.WellType.Injection)
                {
                    q_water[timeLevel] = specifiedFlowRate;
                }
                else
                {
                    // shut-in well
                    //q_oil[time_level] = 0;
                    q_water[timeLevel] = 0;
                    //q_solution_gas[time_level] = 0;
                    //q_water[time_level] = 0;
                }
            }
            else
            {
                if (type == Global.WellType.Production)
                {
                    BHP[timeLevel] = specifiedMinimumBHP;
                    q_oil[timeLevel] = CalculateFlowRate(block.P[timeLevel], BHP[timeLevel], block.Kro[timeLevel], block.viscosityOil[timeLevel], WI, block.Bo[timeLevel]);
                    q_free_gas[timeLevel] = CalculateFlowRate(block.GetPg(timeLevel, timeLevel), BHP[timeLevel], block.Krg[timeLevel], block.viscosityGas[timeLevel], WI, block.Bg[timeLevel]);
                    q_vap_oil[timeLevel] = block.Rvo[timeLevel] * q_free_gas[timeLevel];
                    q_solution_gas[timeLevel] = block.Rso[timeLevel] * q_oil[timeLevel];
                    q_water[timeLevel] = CalculateFlowRate(block.GetPw(timeLevel, timeLevel), BHP[timeLevel], block.Krw[timeLevel], block.viscosityWater[timeLevel], WI, block.Bw[timeLevel]);
                }
                else if (type == Global.WellType.Injection)
                {
                    //To-Do :
                }
                else
                {
                    // shut-in well
                    q_oil[timeLevel] = 0;
                    q_vap_oil[timeLevel] = 0;
                    q_free_gas[timeLevel] = 0;
                    q_solution_gas[timeLevel] = 0;
                    q_water[timeLevel] = 0;
                }
            }
        }

        // checks "for a flow rate controlled well" if the calculated BHP falls below the input minimum value.
        // if so, the well control method changes from rate controlled to specific BHP controlled.
        private bool CheckWellControlChanged(int timeLevel)
        {
            if (control == Global.WellControl.OilRate)
            {
                if (BHP[timeLevel] < specifiedMinimumBHP)
                {
                    control = Global.WellControl.BHP;
                    return true;
                }
            }
            return false;
        }

        // updates the well rates derivatives.
        private void UpdateDerivatives(BaseBlock block)
        {
            double temp;

            if (control == Global.WellControl.OilRate)
            {
                if (type == Global.WellType.Production)
                {
                    // with respect to P
                    dq_oil_dP = 0;
                    temp = CalculatePwf(block.P[2], block.Kro[1], block.viscosityOil[2], block.Bo[2]);
                    dq_free_gas_dP = (CalculateFlowRate(block.GetPg(2, 1), temp, block.Krg[1], block.viscosityGas[2], WI, block.Bg[2]) - q_free_gas[1]) / Global.EPSILON_P;
                    dq_vap_oil_dP = (block.Rvo[2] * CalculateFlowRate(block.P[2], temp, block.Krg[1], block.viscosityGas[2], WI, block.Bg[2]) - block.Rvo[1] * q_free_gas[1]) / Global.EPSILON_P;
                    dq_solution_gas_dP = (block.Rso[2] * CalculateFlowRate(block.P[2], temp, block.Kro[1], block.viscosityOil[2], WI, block.Bo[2]) - block.Rso[1] * q_oil[1]) / Global.EPSILON_P;
                    dq_water_dP = (CalculateFlowRate(block.GetPw(2, 1), temp, block.Krw[1], block.viscosityWater[2], WI, block.Bw[2]) - q_water[1]) / Global.EPSILON_P;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    temp = CalculatePwf(block.P[1], this.data.scal.GetKro(block.Sg[2], block.Sw[1], data.pvt.connateWaterSaturation), block.viscosityOil[1], block.Bo[1]);
                    dq_free_gas_dSg = (CalculateFlowRate(block.GetPg(1, 1), temp, block.Krg[2], block.viscosityGas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.EPSILON_S;
                    dq_vap_oil_dSg = (block.Rvo[1] * CalculateFlowRate(block.P[1], temp, this.data.scal.GetKrg(block.Sg[2]), block.viscosityGas[1], WI, block.Bg[1]) - block.Rvo[1] * q_free_gas[1]) / Global.EPSILON_S;
                    dq_solution_gas_dSg = (block.Rso[1] * /*q_oil[1]*/CalculateFlowRate(block.P[1], temp, this.data.scal.GetKro(block.Sg[2], block.Sw[1], data.pvt.connateWaterSaturation), block.viscosityOil[1], WI, block.Bo[1]) - block.Rso[1] * q_oil[1]) / Global.EPSILON_S;
                    dq_water_dSg = (CalculateFlowRate(block.GetPw(1, 1), temp, block.Krw[1], block.viscosityWater[1], WI, block.Bw[1]) - q_water[1]) / Global.EPSILON_S;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    temp = CalculatePwf(block.P[1], this.data.scal.GetKro(block.Sg[1], block.Sw[2], data.pvt.connateWaterSaturation), block.viscosityOil[1], block.Bo[1]);
                    dq_free_gas_dSw = (CalculateFlowRate(block.GetPg(1, 1), temp, block.Krg[1], block.viscosityGas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.EPSILON_S;
                    dq_vap_oil_dSw = (block.Rvo[1] * CalculateFlowRate(block.P[1], temp, this.data.scal.GetKrg(block.Sg[1]), block.viscosityGas[1], WI, block.Bg[1]) - block.Rvo[1] * q_free_gas[1]) / Global.EPSILON_S;
                    dq_solution_gas_dSw = (block.Rso[1] * CalculateFlowRate(block.P[1], temp, this.data.scal.GetKro(block.Sg[1], block.Sw[2], data.pvt.connateWaterSaturation), block.viscosityOil[2], WI, block.Bo[1]) - block.Rso[1] * q_oil[1]) / Global.EPSILON_S;
                    dq_water_dSw = (CalculateFlowRate(block.GetPw(1, 1), temp, block.Krw[1], block.viscosityWater[1], WI, block.Bw[1]) - q_water[1]) / Global.EPSILON_S;
                }
                else if (type == Global.WellType.Injection)
                {
                    //To-Do : 
                }
                else
                {
                    // shut-in well


                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_vap_oil_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_vap_oil_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_vap_oil_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_water_dSw = 0;
                }

            }
            else if (control == Global.WellControl.GasRate)
            {
                if (type == Global.WellType.Production)
                {
                    //To-Do : 
                }
                else if (type == Global.WellType.Injection)
                {
                    //dq_free_gas_dP = 0;
                    //dq_free_gas_dSg = 0;
                    //dq_free_gas_dSw = 0;


                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_vap_oil_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_vap_oil_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_vap_oil_dSw = 0;
                    dq_water_dSw = 0;
                }
                else
                {
                    // shut-in well


                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_vap_oil_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_vap_oil_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_vap_oil_dSw = 0;
                    dq_water_dSw = 0;
                }
            }
            else if (control == Global.WellControl.WaterRate)
            {
                if (type == Global.WellType.Production)
                {
                    //To-Do : 
                }
                else if (type == Global.WellType.Injection)
                {
                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_vap_oil_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_vap_oil_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_vap_oil_dSw = 0;
                    dq_water_dSw = 0;
                }
                else
                {
                    // shut-in well


                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_vap_oil_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_vap_oil_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_vap_oil_dSw = 0;
                    dq_water_dSw = 0;
                }
            }
            else
            {
                if (type == Global.WellType.Production)
                {
                    temp = BHP[1];
                    // with respect to P
                    dq_oil_dP = (CalculateFlowRate(block.P[2], temp, block.Kro[1], block.viscosityOil[2], WI, block.Bo[2]) - q_oil[1]) / Global.EPSILON_P;
                    dq_free_gas_dP = (CalculateFlowRate(block.GetPg(2, 2), temp, block.Krg[1], block.viscosityGas[2], WI, block.Bg[2]) - q_free_gas[1]) / Global.EPSILON_P;
                    dq_vap_oil_dP = (block.Rvo[2] * CalculateFlowRate(block.GetPg(2, 1), temp, block.Krg[1], block.viscosityGas[2], WI, block.Bg[2]) - block.Rvo[1] * q_free_gas[1]) / Global.EPSILON_P;
                    dq_solution_gas_dP = (block.Rso[2] * CalculateFlowRate(block.P[2], temp, block.Kro[1], block.viscosityOil[2], WI, block.Bo[2]) - block.Rso[1] * q_oil[1]) / Global.EPSILON_P;
                    dq_water_dP = (CalculateFlowRate(block.GetPw(2, 1), temp, block.Krw[1], block.viscosityWater[2], WI, block.Bw[2]) - q_water[1]) / Global.EPSILON_P;
                    // with respect to Sg
                    dq_oil_dSg = (CalculateFlowRate(block.P[1], temp, this.data.scal.GetKro(block.Sg[2], block.Sw[1], data.pvt.connateWaterSaturation), block.viscosityOil[1], WI, block.Bo[1]) - q_oil[1]) / Global.EPSILON_S;
                    dq_free_gas_dSg = (CalculateFlowRate(block.GetPg(1, 1), temp, block.Krg[2], block.viscosityGas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.EPSILON_S;
                    dq_vap_oil_dSg = (block.Rvo[1] * CalculateFlowRate(block.GetPg(1, 2), temp, this.data.scal.GetKrg(block.Sg[2]), block.viscosityGas[1], WI, block.Bg[1]) - block.Rvo[1] * q_free_gas[1]) / Global.EPSILON_S;
                    dq_solution_gas_dSg = (block.Rso[1] * CalculateFlowRate(block.P[1], temp, this.data.scal.GetKro(block.Sg[2], block.Sw[1], data.pvt.connateWaterSaturation), block.viscosityOil[1], WI, block.Bo[1]) - block.Rso[1] * q_oil[1]) / Global.EPSILON_S;
                    dq_water_dSg = (CalculateFlowRate(block.GetPw(1, 1), temp, block.Krw[1], block.viscosityWater[1], WI, block.Bw[1]) - q_water[1]) / Global.EPSILON_S;
                    // with respect to Sw
                    dq_oil_dSw = (CalculateFlowRate(block.P[1], temp, this.data.scal.GetKro(block.Sg[1], block.Sw[2], data.pvt.connateWaterSaturation), block.viscosityOil[1], WI, block.Bo[1]) - q_oil[1]) / Global.EPSILON_S;
                    dq_free_gas_dSw = (CalculateFlowRate(block.GetPg(1, 1), temp, block.Krg[1], block.viscosityGas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.EPSILON_S;
                    dq_vap_oil_dSw = (block.Rvo[1] * CalculateFlowRate(block.GetPg(1, 1), temp, this.data.scal.GetKrg(block.Sg[1]), block.viscosityGas[2], WI, block.Bg[2]) - block.Rvo[1] * q_free_gas[1]) / Global.EPSILON_S;
                    dq_solution_gas_dSw = (block.Rso[1] * CalculateFlowRate(block.P[1], temp, this.data.scal.GetKro(block.Sg[1], block.Sw[2], data.pvt.connateWaterSaturation), block.viscosityOil[2], WI, block.Bo[2]) - block.Rso[1] * q_oil[1]) / Global.EPSILON_S;
                    dq_water_dSw = (CalculateFlowRate(block.GetPw(1, 1), temp, block.Krw[2], block.viscosityWater[1], WI, block.Bw[1]) - q_water[1]) / Global.EPSILON_S;
                }
                else if (type == Global.WellType.Injection)
                {
                    //To-Do :
                }
                else
                {
                    // shut-in well


                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_vap_oil_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_vap_oil_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_vap_oil_dSw = 0;
                    dq_water_dSw = 0;
                }
            }
        }

        /// <summary>
        /// Gets the equivalent well radius.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <param name="grid">The grid.</param>
        /// <returns>The value of the equivalent well radius</returns>
        public static double GetR_Equivalent(BaseBlock block, BaseBlock[] grid)
        {
            double r_equivalnt = 0;

            double distance = 0;
            double numerator = 0, denominator = 0;

            BaseBlock neighbour;
            // neighbour block internal index of "block".
            int index;

            for (int i = 2; i < block.neighborBlocksIndices.Length; i++)
            {
                if (block.neighborBlocksIndices[i] >= 0)
                {
                    neighbour = grid[block.neighborBlocksIndices[i]];
                    index = Array.IndexOf(neighbour.neighborBlocksIndices, block.index);

                    distance = block.deltaXList[i] + neighbour.deltaXList[index];
                    double angle = Math.Atan(0.5 * block.boundaryLengthList[i - 2] / block.deltaXList[i]) * 2;
                    numerator += block.boundaryLengthList[i - 2] / distance * Math.Log(distance) - angle;
                    denominator += (block.boundaryLengthList[i - 2]) / distance;
                }
            }

            r_equivalnt = Math.Exp(numerator / denominator);

            if (grid.Length == 1 || denominator == 0)
            {
                r_equivalnt = 0.2 * block.deltaXList[3] * 2;
            }

            return r_equivalnt;
        }

        /// <summary>
        /// Gets the well productivity index.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <returns></returns>
        private double GetWI(BaseBlock block)
        {
            double K = 0;
            // starting from index = 2, to avoid top and bottom directions permeabilities.
            int counter = 0;
            for (int i = 2; i < block.permeability.Length; i++)
            {
                if (block.permeability[i] > 0)
                {
                    K += block.permeability[i];
                    counter += 1;
                }
                else
                {

                }
            }

            K = K / counter;

            return 2 * Global.Bc * Math.PI * K * block.height / (Math.Log(R_equivalent / radius) + skin);
        }

        /// <summary>
        /// Calculates the bottom hole flowing pressure.
        /// </summary>
        /// <param name="pressure">The pressure.</param>
        /// <param name="Kr">The kr.</param>
        /// <param name="viscosity">The viscosity.</param>
        /// <param name="FVF">The FVF.</param>
        /// 
        /// <returns>The value of PWF.</returns>
        public double CalculatePwf(double pressure, double Kr, double viscosity, double FVF)
        {
            double mobility = Kr / (viscosity * FVF);
            return pressure - (q_oil[0] / (mobility * WI));
        }

        /// <summary>
        /// Calculates the flow rate of a certain phase.
        /// </summary>
        /// <param name="P">The p.</param>
        /// <param name="Pwf">The PWF.</param>
        /// <param name="Kr">The kr.</param>
        /// <param name="viscosity">The viscosity.</param>
        /// <param name="WI">The wi.</param>
        /// <param name="FVF">The FVF.</param>
        /// <returns>The calculated flow rate</returns>
        public static double CalculateFlowRate(double P, double Pwf, double Kr, double viscosity, double WI, double FVF)
        {
            return (P - Pwf) * (Kr / (viscosity * FVF)) * WI;
            //return 0;
        }

        /// <summary>
        /// Gets the well mobility of a certain phase at a specified time level.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <param name="phase">The phase.</param>
        /// <param name="timeLevel">The time level.</param>
        /// <returns></returns>
        /// <seealso cref="Global.STEPS_MEMORY"/>
        public double GetWellMobility(SimulationData data, Global.Phase phase, int timeLevel)
        {
            BaseBlock block = data.grid[this.index];

            if (phase == Global.Phase.Oil)
            {
                return block.Kro[timeLevel] / (block.viscosityOil[timeLevel] * block.Bo[timeLevel]);
            }
            else if (phase == Global.Phase.Gas)
            {
                return block.Krg[timeLevel] / (block.viscosityGas[timeLevel] * block.Bg[timeLevel]);
            }
            else
            {
                return block.Krw[timeLevel] / (block.viscosityWater[timeLevel] * block.Bw[timeLevel]);
            }
        }

        /// <summary>
        /// Gets the well mobility according to the values of Kr, viscosity and formation volume factor input to the method.
        /// </summary>
        /// <remarks>
        /// This method differs from <see cref="GetWellMobility(SimulationData, Global.Phase, int)"/> in that it specifies the PVT and SCAL values instead of infering them automatically from the phase.
        /// </remarks>
        /// <param name="Kr">The kr.</param>
        /// <param name="viscosity">The viscosity.</param>
        /// <param name="FVF">The FVF.</param>
        /// <returns></returns>
        public static double GetWellMobility(double Kr, double viscosity, double FVF)
        {
            return Kr / (viscosity * FVF);
        }


        public double GetProducingGOR()
        {
            return (q_solution_gas[0] + q_free_gas[0]) * Global.a / q_oil[0] / 1000;
        }

        public double GetTotalOilProduction()
        {
            return q_oil[0] + q_vap_oil[0];
        }

        public double GetProductionOil()
        {
            return q_oil[0] + q_vap_oil[0];
        }
    }
}
