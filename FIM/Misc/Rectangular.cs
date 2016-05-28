using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FIM.Misc
{
    class Rectangular
    {
        public static int xyzToNatural(int x, int y, int z, int i, int j, int k)
        {
            return i + j * x + k * x * y;
        }
    }
}
