using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace ConsoleApp1
{
    class Program
    {
        static void Main(string[] args)
        {
            string filename = args[0];

            Dictionary<string, HashSet<string>> subsets = new Dictionary<string, HashSet<string>>();
            using (StreamReader input = new StreamReader(filename))
            {
                while (!input.EndOfStream)
                {
                    string subset = input.ReadLine();
                    string line = input.ReadLine();
                    string[] fields = line.Split(',');
                    subsets.Add(subset, new HashSet<string>(fields));
                }
            }

            Console.WriteLine(',' + string.Join(",", subsets.Keys));
            foreach (string isoform in subsets["overall"])
            {
                StringBuilder sb = new StringBuilder(isoform + ',');
                foreach (KeyValuePair<string, HashSet<string>> kvp in subsets)
                    sb.Append((kvp.Value.Contains(isoform) ? 1 : 0).ToString() + ',');
                sb.Remove(sb.Length - 1, 1);
                Console.WriteLine(sb.ToString());
            }
        }
    }
}
