using FIM.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Well
{
    /// <summary>
    /// This class contains helper methods for handling well calculations.
    /// </summary>
    public class WellData
    {
        public int index;
        public double radius;
        public double skin;

        public Global.WellType type;
        public Global.WellControl control;
        public double specified_minimum_BHP;
        public double specified_flow_rate;
        public double R_equivalent;
        public double WI;
        public double[] BHP;
        public double[] q_oil, q_free_gas, q_solution_gas, q_water;

        // derivatives
        public double dq_oil_dP, dq_oil_dSg, dq_oil_dSw;
        public double dq_free_gas_dP, dq_free_gas_dSg, dq_free_gas_dSw;
        public double dq_solution_gas_dP, dq_solution_gas_dSg, dq_solution_gas_dSw;
        public double dq_water_dP, dq_water_dSg, dq_water_dSw;

        public WellData(SimulationData data, int index, Global.WellType type, Global.WellControl control, double radius, double skin, double specified_minimum_BHP = 0, double specified_flow_rate = 0)
        {
            this.index = index;
            data.grid[index].type = Global.BlockType.Well_Block;
            this.radius = radius;
            this.skin = skin;
            this.R_equivalent = getR_equivalent(data.grid[index], data.grid);
            this.WI = getWI(data.grid[index]);
            this.specified_minimum_BHP = specified_minimum_BHP;
            
            this.specified_flow_rate = specified_flow_rate;
            if (type == Global.WellType.Injection)
            {
                this.specified_flow_rate = -1 * specified_flow_rate;
            }

            this.BHP = new double[Global.steps_memory];
            this.q_oil = new double[Global.steps_memory];
            this.q_free_gas = new double[Global.steps_memory];
            this.q_solution_gas = new double[Global.steps_memory];
            this.q_water = new double[Global.steps_memory];

            this.type = type;
            this.control = control;
        }

        public void update(int time_level, BaseBlock block)
        {
            updateFlowRates(time_level, block);
            if (checkWellControlChanged(time_level))
            {
                updateFlowRates(time_level, block);
            }
            //checkWellControlChanged();
            updateDerivatives(block);
        }

        private void updateFlowRates(int time_level, BaseBlock block)
        {
            if (control == Global.WellControl.OilRate)
            {
                if (type == Global.WellType.Production)
                {
                    q_oil[time_level] = specified_flow_rate;
                    BHP[time_level] = calculatePwf(block, block.P[time_level], block.Kro[time_level], block.viscosity_oil[time_level], block.Bo[time_level]);
                    q_free_gas[time_level] = calculateFlow_Rate(block.P[time_level], BHP[time_level], block.Krg[time_level], block.viscosity_gas[time_level], WI, block.Bg[time_level]);
                    q_solution_gas[time_level] = block.Rso[time_level] * q_oil[time_level];
                    q_water[time_level] = calculateFlow_Rate(block.P[time_level], BHP[time_level], block.Krw[time_level], block.viscosity_water[time_level], WI, block.Bw[time_level]);
                }
                else if (type == Global.WellType.Injection)
                {
                    //To-Do : 
                }
                else
                {
                    // shut-in well
                    q_oil[time_level] = 0;
                    q_free_gas[time_level] = 0;
                    q_solution_gas[time_level] = 0;
                    q_water[time_level] = 0;
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
                    q_free_gas[time_level] = specified_flow_rate;
                }
                else
                {
                    // shut-in well
                    //q_oil[time_level] = 0;
                    q_free_gas[time_level] = 0;
                    //q_solution_gas[time_level] = 0;
                    //q_water[time_level] = 0;
                }
            }
            else
            {
                if (type == Global.WellType.Production)
                {
                    BHP[time_level] = specified_minimum_BHP;
                    q_oil[time_level] = calculateFlow_Rate(block.P[time_level], BHP[time_level], block.Kro[time_level], block.viscosity_oil[time_level], WI, block.Bo[time_level]);
                    q_free_gas[time_level] = calculateFlow_Rate(block.P[time_level], BHP[time_level], block.Krg[time_level], block.viscosity_gas[time_level], WI, block.Bg[time_level]);
                    q_solution_gas[time_level] = block.Rso[time_level] * q_oil[time_level];
                    q_water[time_level] = calculateFlow_Rate(block.P[time_level], BHP[time_level], block.Krw[time_level], block.viscosity_water[time_level], WI, block.Bw[time_level]);
                }
                else if (type == Global.WellType.Injection)
                {
                    //To-Do :
                }
                else
                {
                    // shut-in well
                    q_oil[time_level] = 0;
                    q_free_gas[time_level] = 0;
                    q_solution_gas[time_level] = 0;
                    q_water[time_level] = 0;
                }
            }
        }

        private bool checkWellControlChanged(int time_level)
        {
            if (control == Global.WellControl.OilRate)
            {
                if (BHP[time_level] < specified_minimum_BHP)
                {
                    control = Global.WellControl.BHP;
                    return true;
                }
            }
            return false;
        }

        private void updateDerivatives(BaseBlock block)
        {
            double temp;

            if (control == Global.WellControl.OilRate)
            {
                if (type == Global.WellType.Production)
                {
                    // with respect to P
                    dq_oil_dP = 0;
                    temp = calculatePwf(block, block.P[2], block.Kro[1], block.viscosity_oil[2], block.Bo[2]);
                    dq_free_gas_dP = (calculateFlow_Rate(block.P[2], temp, block.Krg[1], block.viscosity_gas[2], WI, block.Bg[2]) - q_free_gas[1]) / Global.epsilon;
                    dq_solution_gas_dP = (block.Rso[2] * calculateFlow_Rate(block.P[2], temp, block.Kro[1], block.viscosity_oil[2], WI, block.Bo[2]) - block.Rso[1] * q_oil[1]) / Global.epsilon;
                    dq_water_dP = (calculateFlow_Rate(block.P[2], temp, block.Krw[1], block.viscosity_water[2], WI, block.Bw[2]) - q_water[1]) / Global.epsilon;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    temp = calculatePwf(block, block.P[1], block.Kro[2], block.viscosity_oil[1], block.Bo[1]);
                    dq_free_gas_dSg = (calculateFlow_Rate(block.P[1], temp, block.Krg[2], block.viscosity_gas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.epsilon;
                    dq_solution_gas_dSg = (block.Rso[1] * q_oil[1] - block.Rso[1] * q_oil[1]) / Global.epsilon;
                    dq_water_dSg = (calculateFlow_Rate(block.P[1], temp, block.Krw[1], block.viscosity_water[1], WI, block.Bw[1]) - q_water[1]) / Global.epsilon;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    temp = calculatePwf(block, block.P[1], block.Kro[1], block.viscosity_oil[1], block.Bo[1]);
                    dq_free_gas_dSw = (calculateFlow_Rate(block.P[1], temp, block.Krg[1], block.viscosity_gas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.epsilon;
                    dq_solution_gas_dSw = (block.Rso[1] * q_oil[1] - block.Rso[1] * q_oil[1]) / Global.epsilon;
                    dq_water_dSw = (calculateFlow_Rate(block.P[1], temp, block.Krw[1], block.viscosity_water[1], WI, block.Bw[1]) - q_water[1]) / Global.epsilon;
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
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
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
                    dq_free_gas_dP = 0;
                    dq_free_gas_dSg = 0;
                    dq_free_gas_dSw = 0;
                }
                else
                {
                    // shut-in well


                    // with respect to P
                    dq_oil_dP = 0;
                    dq_free_gas_dP = 0;
                    dq_solution_gas_dP = 0;
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
                    dq_water_dSw = 0;
                }
            }
            else
            {
                if (type == Global.WellType.Production)
                {
                    temp = BHP[1];
                    // with respect to P
                    dq_oil_dP = (calculateFlow_Rate(block.P[2], temp, block.Kro[1], block.viscosity_oil[2], WI, block.Bo[2]) - q_oil[1]) / Global.epsilon;
                    dq_free_gas_dP = (calculateFlow_Rate(block.P[2], temp, block.Krg[1], block.viscosity_gas[2], WI, block.Bg[2]) - q_free_gas[1]) / Global.epsilon;
                    dq_solution_gas_dP = (block.Rso[2] * calculateFlow_Rate(block.P[2], temp, block.Kro[1], block.viscosity_oil[2], WI, block.Bo[2]) - block.Rso[1] * q_oil[1]) / Global.epsilon;
                    dq_water_dP = (calculateFlow_Rate(block.P[2], temp, block.Krw[1], block.viscosity_water[2], WI, block.Bw[2]) - q_water[1]) / Global.epsilon;
                    // with respect to Sg
                    dq_oil_dSg = (calculateFlow_Rate(block.P[1], temp, block.Kro[2], block.viscosity_oil[1], WI, block.Bo[1]) - q_oil[1]) / Global.epsilon;
                    dq_free_gas_dSg = (calculateFlow_Rate(block.P[1], temp, block.Krg[2], block.viscosity_gas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.epsilon;
                    dq_solution_gas_dSg = (block.Rso[1] * q_oil[1] - block.Rso[1] * q_oil[1]) / Global.epsilon;
                    dq_water_dSg = (calculateFlow_Rate(block.P[1], temp, block.Krw[1], block.viscosity_water[1], WI, block.Bw[1]) - q_water[1]) / Global.epsilon;
                    // with respect to Sw
                    dq_oil_dSw = (calculateFlow_Rate(block.P[1], temp, block.Kro[1], block.viscosity_oil[1], WI, block.Bo[1]) - q_oil[1]) / Global.epsilon;
                    dq_free_gas_dSw = (calculateFlow_Rate(block.P[1], temp, block.Krg[2], block.viscosity_gas[1], WI, block.Bg[1]) - q_free_gas[1]) / Global.epsilon;
                    dq_solution_gas_dSw = (block.Rso[1] * q_oil[1] - block.Rso[1] * q_oil[1]) / Global.epsilon;
                    dq_water_dSw = (calculateFlow_Rate(block.P[1], temp, block.Krw[1], block.viscosity_water[1], WI, block.Bw[1]) - q_water[1]) / Global.epsilon;
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
                    dq_water_dP = 0;
                    // with respect to Sg
                    dq_oil_dSg = 0;
                    dq_free_gas_dSg = 0;
                    dq_solution_gas_dSg = 0;
                    dq_water_dSg = 0;
                    // with respect to Sw
                    dq_oil_dSw = 0;
                    dq_free_gas_dSw = 0;
                    dq_solution_gas_dSw = 0;
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
        public static double getR_equivalent(BaseBlock block, BaseBlock[] grid)
        {
            double r_equivalnt = 0;

            double distance = 0;
            double numerator = 0, denominator = 0;

            BaseBlock neighbour;
            // neighbour block internal index of "block".
            int index;

            for (int i = 2; i < block.neighbour_blocks_indices.Length; i++)
            {
                if (block.neighbour_blocks_indices[i] >= 0)
                {
                    neighbour = grid[block.neighbour_blocks_indices[i]];
                    index = Array.IndexOf(neighbour.neighbour_blocks_indices, block.index);

                    distance = (block.delta_x_list[i] + neighbour.delta_x_list[index]) / 2;
                    numerator += block.boundary_length_list[i] / distance * Math.Log(distance) - 0.5 * Math.PI;
                    denominator += (block.boundary_length_list[i]) / distance;
                }
            }

            r_equivalnt = Math.Exp(numerator / denominator);

            if (grid.Length == 1)
            {
                r_equivalnt = 0.2 * block.delta_x_list[3];
            }

            return r_equivalnt;
        }

        /// <summary>
        /// Gets the well productivity index.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <returns></returns>
        private double getWI(BaseBlock block)
        {
            double K = 0;
            // starting from index = 2, to avoid top and bottom directions permeabilities.
            int counter = 0;
            for (int i = 2; i < block.permeability_list.Length; i++)
            {
                if (block.permeability_list[i] > 0)
                {
                    K += block.permeability_list[i];
                    counter += 1;
                }
                else
                {
                    
                }
            }

            K = K / counter;

            return 2 * Global.Bc * Math.PI * K * block.h / (Math.Log(R_equivalent / radius) + skin);
        }

        /// <summary>
        /// Calculates the bottom hole flowing pressure.
        /// </summary>
        /// <param name="block">The block.</param>
        /// <param name="pressure">The pressure.</param>
        /// <param name="Kr">The kr.</param>
        /// <param name="viscosity">The viscosity.</param>
        /// <param name="FVF">The FVF.</param>
        /// <returns>The value of PWF.</returns>
        private double calculatePwf(BaseBlock block, double pressure, double Kr, double viscosity, double FVF)
        {
            double mobility = Kr / (viscosity * FVF);
            return pressure - ( q_oil[0] / (mobility * WI) );
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
        public static double calculateFlow_Rate(double P, double Pwf, double Kr, double viscosity, double WI, double FVF)
        {
            return (P - Pwf) * (Kr / (viscosity * FVF)) * WI;
            //return 0;
        }
    }
}
