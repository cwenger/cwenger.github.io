using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace MaxQuantAnalyzer3
{
    class Program
    {
        static readonly Regex EXPERIMENT_MATCH = new Regex(@"(.+)_(.+)_(.+)");

        static void Main(string[] args)
        {
            string peptides_filename = args[0];
            bool exclude_special_proteins = true;

            // read in peptides.txt and make a list of all the peptides containg the sequence and the number of times it's found in each experiment, keyed by peptide ID
            Dictionary<string, Counts> subsets = new Dictionary<string, Counts>();
            Counts overall = new Counts();
            using (StreamReader peptides_file = new StreamReader(peptides_filename))
            {
                string header = peptides_file.ReadLine();
                string[] header_fields = header.Split('\t');
                int rev_index = Array.IndexOf(header_fields, "Reverse");
                int con_index = Array.IndexOf(header_fields, "Potential contaminant");
                Dictionary<int, string> experiment_indexes = new Dictionary<int, string>();
                for (int h = 0; h < header_fields.Length; h++)
                {
                    // make a list of each possible subset of runs
                    // also create a list of all experiments keyed by their column index
                    if (header_fields[h].StartsWith("Experiment"))
                    {
                        string experiment = header_fields[h].Substring("Experiment ".Length);
                        Match match = EXPERIMENT_MATCH.Match(experiment);
                        if (!subsets.ContainsKey(match.Groups[1].Value))
                            subsets.Add(match.Groups[1].Value, new Counts());
                        if (!subsets.ContainsKey(match.Groups[2].Value))
                            subsets.Add(match.Groups[2].Value, new Counts());
                        if (!subsets.ContainsKey(match.Groups[3].Value))
                            subsets.Add(match.Groups[3].Value, new Counts());
                        if (!subsets.ContainsKey(match.Groups[1].Value + '_' + match.Groups[2].Value))
                            subsets.Add(match.Groups[1].Value + '_' + match.Groups[2].Value, new Counts());
                        if (!subsets.ContainsKey(match.Groups[2].Value + '_' + match.Groups[3].Value))
                            subsets.Add(match.Groups[2].Value + '_' + match.Groups[3].Value, new Counts());
                        if (!subsets.ContainsKey(match.Groups[1].Value + '_' + match.Groups[3].Value))
                            subsets.Add(match.Groups[1].Value + '_' + match.Groups[3].Value, new Counts());
                        subsets.Add(experiment, new Counts());
                        experiment_indexes.Add(h, experiment);
                    }
                }

                while (!peptides_file.EndOfStream)
                {
                    string line = peptides_file.ReadLine();
                    string[] fields = line.Split('\t');

                    if (exclude_special_proteins && (!string.IsNullOrWhiteSpace(fields[rev_index]) || !string.IsNullOrWhiteSpace(fields[con_index])))
                        continue;

                    Dictionary<string, int> psms_per_exp = new Dictionary<string, int>();
                    foreach (KeyValuePair<int, string> kvp in experiment_indexes)
                    {
                        int.TryParse(fields[kvp.Key], out int psms);
                        psms_per_exp.Add(kvp.Value, psms);
                    }

                    foreach (string subset in subsets.Keys)
                    {
                        string[] subset_components = subset.Split('_');

                        int psms = 0;
                        foreach (KeyValuePair<string, int> kvp in psms_per_exp)
                            if (Array.TrueForAll(subset_components, x => kvp.Key.Contains(x)))
                                psms += kvp.Value;
                        subsets[subset].PSMs += psms;
                        subsets[subset].UniquePeptides += psms > 0 ? 1 : 0;
                    }

                    int overall_psms = 0;
                    foreach (KeyValuePair<string, int> kvp in psms_per_exp)
                        overall_psms += kvp.Value;
                    overall.PSMs += overall_psms;
                    overall.UniquePeptides++;
                }
            }

            Console.WriteLine("Subset\tPSMs\tUnique Peptides");
            foreach (KeyValuePair<string, Counts> kvp in subsets)
                Console.WriteLine("{0}\t{1}\t{2}", kvp.Key, kvp.Value.PSMs, kvp.Value.UniquePeptides);
            Console.WriteLine("{0}\t{1}\t{2}", "overall", overall.PSMs, overall.UniquePeptides);
            Console.ReadKey(true);
        }

        public class Counts
        {
            public int PSMs { get; set; }
            public int UniquePeptides { get; set; }
        }
    }
}
