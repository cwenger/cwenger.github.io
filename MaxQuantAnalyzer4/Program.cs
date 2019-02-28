using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace MaxQuantAnalyzer4
{
    class Program
    {
        static readonly Regex EXPERIMENT_MATCH = new Regex(@"(.+)_(.+)_(.+)");

        static void Main(string[] args)
        {
            string proteins_filename = args[0];
            bool exclude_special_proteins = true;

            Dictionary<string, Counts> subsets = new Dictionary<string, Counts>();
            Counts overall = new Counts();
            using (StreamReader peptides_file = new StreamReader(proteins_filename))
            {
                string header = peptides_file.ReadLine();
                string[] header_fields = header.Split('\t');
                int site_only_index = Array.IndexOf(header_fields, "Only identified by site");
                int rev_index = Array.IndexOf(header_fields, "Reverse");
                int con_index = Array.IndexOf(header_fields, "Potential contaminant");
                Dictionary<int, string> peptide_indexes = new Dictionary<int, string>();
                Dictionary<int, string> razor_plus_unique_peptide_indexes = new Dictionary<int, string>();
                Dictionary<int, string> unique_peptide_indexes = new Dictionary<int, string>();
                for (int h = 0; h < header_fields.Length; h++)
                {
                    // make a list of each possible subset of runs
                    // also create a list of all experiments keyed by their column index
                    if (header_fields[h].StartsWith("Peptides "))
                    {
                        string experiment = header_fields[h].Substring("Peptides ".Length);
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
                        peptide_indexes.Add(h, experiment);
                    }
                    else if (header_fields[h].StartsWith("Razor + unique peptides "))
                    {
                        string experiment = header_fields[h].Substring("Razor + unique peptides ".Length);
                        razor_plus_unique_peptide_indexes.Add(h, experiment);
                    }
                    else if (header_fields[h].StartsWith("Unique peptides "))
                    {
                        string experiment = header_fields[h].Substring("Unique peptides ".Length);
                        unique_peptide_indexes.Add(h, experiment);
                    }
                }

                while (!peptides_file.EndOfStream)
                {
                    string line = peptides_file.ReadLine();
                    string[] fields = line.Split('\t');

                    if (exclude_special_proteins && (fields[0].StartsWith("REV") || fields[0].StartsWith("CON") || !string.IsNullOrWhiteSpace(fields[site_only_index]) || !string.IsNullOrWhiteSpace(fields[rev_index]) || !string.IsNullOrWhiteSpace(fields[con_index])))
                        continue;

                    Dictionary<string, int> peptides = new Dictionary<string, int>();
                    foreach (KeyValuePair<int, string> kvp in peptide_indexes)
                    {
                        int.TryParse(fields[kvp.Key], out int peptide_count);
                        peptides.Add(kvp.Value, peptide_count);
                    }
                    Dictionary<string, int> razor_plus_unique_peptides = new Dictionary<string, int>();
                    foreach (KeyValuePair<int, string> kvp in razor_plus_unique_peptide_indexes)
                    {
                        int.TryParse(fields[kvp.Key], out int razor_plus_unique_peptide_count);
                        razor_plus_unique_peptides.Add(kvp.Value, razor_plus_unique_peptide_count);
                    }
                    Dictionary<string, int> unique_peptides = new Dictionary<string, int>();
                    foreach (KeyValuePair<int, string> kvp in unique_peptide_indexes)
                    {
                        int.TryParse(fields[kvp.Key], out int unique_peptide_count);
                        unique_peptides.Add(kvp.Value, unique_peptide_count);
                    }

                    foreach (string subset in subsets.Keys)
                    {
                        string[] subset_components = subset.Split('_');

                        int peptide_count = 0;
                        foreach (KeyValuePair<string, int> kvp in peptides)
                            if (Array.TrueForAll(subset_components, x => kvp.Key.Contains(x)))
                                peptide_count += kvp.Value;
                        int razor_plus_unique_peptide_count = 0;
                        foreach (KeyValuePair<string, int> kvp in razor_plus_unique_peptides)
                            if (Array.TrueForAll(subset_components, x => kvp.Key.Contains(x)))
                                razor_plus_unique_peptide_count += kvp.Value;
                        int unique_peptide_count = 0;
                        foreach (KeyValuePair<string, int> kvp in unique_peptides)
                            if (Array.TrueForAll(subset_components, x => kvp.Key.Contains(x)))
                                unique_peptide_count += kvp.Value;
                        subsets[subset].Peptides += peptide_count;
                        subsets[subset].RazorPlusUniquePeptides += razor_plus_unique_peptide_count;
                        subsets[subset].UniquePeptides += unique_peptide_count;
                        subsets[subset].ProteinGroups += peptide_count > 0 ? 1 : 0;
                    }

                    int overall_peptide_count = 0;
                    foreach (KeyValuePair<string, int> kvp in peptides)
                        overall_peptide_count += kvp.Value;
                    int overall_razor_plus_unique_peptide_count = 0;
                    foreach (KeyValuePair<string, int> kvp in razor_plus_unique_peptides)
                        overall_razor_plus_unique_peptide_count += kvp.Value;
                    int overall_unique_peptide_count = 0;
                    foreach (KeyValuePair<string, int> kvp in unique_peptides)
                        overall_unique_peptide_count += kvp.Value;
                    overall.Peptides += overall_peptide_count;
                    overall.RazorPlusUniquePeptides += overall_razor_plus_unique_peptide_count;
                    overall.UniquePeptides += overall_unique_peptide_count;
                    overall.ProteinGroups++;
                }
            }

            Console.WriteLine("Subset\tPeptides\tRazor + Unique Peptides\tUnique Peptides\tProtein Groups");
            foreach (KeyValuePair<string, Counts> kvp in subsets)
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kvp.Key, kvp.Value.Peptides, kvp.Value.RazorPlusUniquePeptides, kvp.Value.UniquePeptides, kvp.Value.ProteinGroups);
            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", "overall", overall.Peptides, overall.RazorPlusUniquePeptides, overall.UniquePeptides, overall.ProteinGroups);
            Console.ReadKey(true);
        }

        public class Counts
        {
            public int Peptides { get; set; }
            public int RazorPlusUniquePeptides { get; set; }
            public int UniquePeptides { get; set; }
            public int ProteinGroups { get; set; }
        }
    }
}
