using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CrashingCost
{
    class Diagram
    {
        public List<int> critical_path_indices;

        public Activity[] activity_list;

        public double total_project_time;

        public Diagram(Activity[] activity_list)
        {
            this.activity_list = activity_list;

            setSuccessors();
        }

        private void setSuccessors()
        {
            Activity activity;

            for (int i = 0; i < activity_list.Length; i++)
            {
                for (int j = 0; j < activity_list[i].predecessors_list.Count; j++)
                {
                    activity = activity_list[activity_list[i].predecessors_list[j]];
                    if (!activity.successors_list.Contains(i))
                    {
                        activity.successors_list.Add(i);
                    }
                }
            }
        }

        public void calculateES_EF()
        {
            Activity activity;

            for (int i = 0; i < this.activity_list.Length; i++)
            {
                activity = activity_list[i];
                activity.setES_EF(activity_list);
            }

            this.total_project_time = activity_list.Select(x => x.EF).Max();
        }

        public void calculateLS_LF()
        {
            Activity activity;

            for (int i = 0; i < this.activity_list.Length; i++)
            {
                activity = activity_list[i];
                activity.setLS_LF(activity_list, total_project_time);
            }
        }

        public void calculateCriticalPath()
        {
            List<int> critical_activities = new List<int>();

            for (int i = 0; i < activity_list.Length; i++)
            {
                if (activity_list[i].ES == activity_list[i].LS)
                {
                    critical_activities.Add(i);
                }
            }

            this.critical_path_indices = critical_activities;
            
        }
    }
}
