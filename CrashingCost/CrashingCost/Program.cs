using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrashingCost
{
    class Program
    {
        static void Main(string[] args)
        {
            Activity A = new Activity("A", 14, 1);
            Activity B = new Activity("B", 3, 1);
            Activity C = new Activity("C", 7, 1);
            Activity D = new Activity("D", 4, 1);
            Activity E = new Activity("E", 3, 1);
            Activity F = new Activity("F", 10, 1);

            Activity[] activities = new Activity[] { A, B, C, D, E, F };

            C.setPredecessors(activities, A);
            E.setPredecessors(activities, B);
            D.setPredecessors(activities, C);
            F.setPredecessors(activities, D, E);

            Diagram diagram = new Diagram(activities);
            diagram.calculateES_EF();
            diagram.calculateLS_LF();
            diagram.calculateCriticalPath();

            Console.ReadKey();
        }
    }
}
