using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace FIM.Parser
{
    public static class KeyWordMatcher
    {
        /// <summary>
        /// Determines whether the line contains any of the specified key words.
        /// </summary>
        /// <param name="line">The line.</param>
        /// <param name="KeyWords">The key words.</param>
        /// <returns></returns>
        public static bool ContainsAnyOf(this string line, string[] KeyWords)
        {
            for (int i = 0; i < KeyWords.Length; i++)
            {
                //if (line.Contains(KeyWords[i]))
                //{
                //    return true;
                //}

                if (Regex.Match(line, "\\b" + KeyWords[i] + "\\b").Success)
                {
                    return true;
                }
            }

            return false;
        }
    }
}
