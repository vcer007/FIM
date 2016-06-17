using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrashingCost
{
    class Activity
    {
        public bool ES_EF = false;
        public bool LS_LF = false;

        public string activity_name;
        public double duration;
        public double cost_slope;


        public List<int> predecessors_list;
        public List<int> successors_list;

        public double ES, EF, LS, LF, Total_Float;

        public double median, variance;

        public Activity(string name, double duration, double cost_slope)
        {
            this.activity_name = name;
            this.duration = duration;
            this.cost_slope = cost_slope;

            predecessors_list = new List<int>();
            successors_list = new List<int>();
        }

        public void setES_EF(Activity[] activity_list)
        {
            this.ES_EF = true;

            if (this.predecessors_list.Count == 0)
            {
                this.ES = 0;
                this.EF = this.duration;
            }
            else
            {
                for (int i = 0; i < this.predecessors_list.Count; i++)
                {
                    if (!activity_list[predecessors_list[i]].ES_EF)
                    {
                        activity_list[predecessors_list[i]].setES_EF(activity_list);
                    }
                }

                List<double> temp = new List<double>();
                for (int i = 0; i < predecessors_list.Count; i++)
                {
                    temp.Add(activity_list[predecessors_list[i]].EF);
                }
                this.ES = temp.Max();
                this.EF = this.ES + this.duration;
            }
        }

        public void setLS_LF(Activity[] activity_list, double end_time)
        {
            this.LS_LF = true;

            if (this.successors_list.Count == 0)
            {
                this.LF = end_time;
                this.LS = this.LF - this.duration;
            }
            else
            {
                for (int i = 0; i < this.successors_list.Count; i++)
                {
                    if (!activity_list[successors_list[i]].LS_LF)
                    {
                        activity_list[successors_list[i]].setLS_LF(activity_list, end_time);
                    }
                }

                this.LF = this.successors_list.Select(x => activity_list[x].LS).Min();
                this.LS = this.LF - this.duration;
            }
        }

        public void setPredecessors(Activity[] activities_list, params Activity[] activity)
        {
            for (int i = 0; i < activity.Length; i++)
            {
                this.predecessors_list.Add(Array.IndexOf(activities_list, activity[i]));
            }
        }
    }
}
