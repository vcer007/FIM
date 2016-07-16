using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Rock
{
    public class Porosity
    {
        //Rock compressibility factor
        double Cf;
        double porosity_ref;
        double pressure_ref;

        //Constructor that takes a Cf value as an argument
        public Porosity(double Cf, double porosity_ref, double pressure_ref)
        {
            this.Cf = Cf;
            this.porosity_ref = porosity_ref;
            this.pressure_ref = pressure_ref;
        }

        //The method used to generate the new pore volume
        public double getPorosity(double pressure)
        {
            return porosity_ref * (1 + Cf * (pressure - pressure_ref));
        }

    }
}
