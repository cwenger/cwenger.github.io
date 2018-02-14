using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace MaxQuantAnalyzer2
{
    class Program
    {
        static readonly Regex EXPERIMENT_MATCH = new Regex(@"(.+)_(.+)_(.+)");

        enum AnalysisType
        {
            Independent,
            ConserveGroups
        }

        const AnalysisType ANALYSIS_TYPE = AnalysisType.ConserveGroups;
        const bool REMOVE_ISOFORMS_WITHOUT_UNIQUE_PEPTIDES = false;
        const bool REMOVE_EACH_GROUP_SEQUENTIALLY_TEST = false;

        static void Main(string[] args)
        {
            string peptides_filename = args[0];
            bool exclude_nomsms_peptides;
            if (args.Length >= 2)
                bool.TryParse(args[1], out exclude_nomsms_peptides);
            else
                exclude_nomsms_peptides = false;

            List<Peptide> all_peptides = new List<Peptide>();
            Dictionary<string, List<Peptide>> isoforms = new Dictionary<string, List<Peptide>>();
            HashSet<string> subsets = new HashSet<string>();
            using (StreamReader peptides_file = new StreamReader(peptides_filename))
            {
                string header = peptides_file.ReadLine();
                string[] header_fields = header.Split('\t');
                int id_index = Array.IndexOf(header_fields, "id");
                int best_msms_index = Array.IndexOf(header_fields, "Best MS/MS");
                int score_index = Array.IndexOf(header_fields, "Score");
                int isoforms_index = Array.IndexOf(header_fields, "isoforms(+)");
                Dictionary<int, string> experiment_indexes = new Dictionary<int, string>();
                for (int h = 0; h < header_fields.Length; h++)
                {
                    // make a list of each possible subset of runs
                    // also create a list of all experiments keyed by their column index
                    if (header_fields[h].StartsWith("Experiment"))
                    {
                        string experiment = header_fields[h].Substring("Experiment ".Length);
                        Match match = EXPERIMENT_MATCH.Match(experiment);
                        subsets.Add(match.Groups[1].Value);  // only care about cell lines
                        //subsets.Add(match.Groups[2].Value);
                        //subsets.Add(match.Groups[3].Value);
                        //subsets.Add(match.Groups[1].Value + '_' + match.Groups[2].Value);
                        //subsets.Add(match.Groups[2].Value + '_' + match.Groups[3].Value);
                        //subsets.Add(match.Groups[1].Value + '_' + match.Groups[3].Value);
                        //subsets.Add(experiment);
                        experiment_indexes.Add(h, experiment);
                    }
                }

                while (!peptides_file.EndOfStream)
                {
                    string line = peptides_file.ReadLine();
                    string[] fields = line.Split('\t');

                    if (string.IsNullOrWhiteSpace(fields[isoforms_index]))
                        continue;

                    if (exclude_nomsms_peptides && string.IsNullOrWhiteSpace(fields[best_msms_index]))
                        continue;

                    HashSet<string> experiments = new HashSet<string>();
                    foreach (KeyValuePair<int, string> kvp in experiment_indexes)
                    {
                        int psms;
                        int.TryParse(fields[kvp.Key], out psms);
                        if (psms > 0)
                            experiments.Add(kvp.Value);
                    }

                    Peptide peptide = new Peptide(int.Parse(fields[id_index]), double.Parse(fields[score_index]), experiments);
                    all_peptides.Add(peptide);
                    foreach (string isoform in fields[isoforms_index].Split(';'))
                    {
                        List<Peptide> peptides;
                        if (!isoforms.TryGetValue(isoform, out peptides))
                        {
                            peptides = new List<Peptide>();
                            peptides.Add(peptide);
                            isoforms.Add(isoform, peptides);
                        }
                        else
                            peptides.Add(peptide);
                    }
                }
            }

            List<KeyValuePair<List<string>, List<Peptide>>> overall_minimal_isoforms = FindMinimalProteinGroups(isoforms);
            Console.Error.WriteLine("overall" + '\t' + overall_minimal_isoforms.Count.ToString());

            HashSet<int> peptide_ids_from_parsimony_isoforms = new HashSet<int>(overall_minimal_isoforms.SelectMany(x => x.Value.Select(y => y.Id)));
            if (!new HashSet<int>(all_peptides.Select(x => x.Id)).SetEquals(peptide_ids_from_parsimony_isoforms))
                throw new Exception();

            Console.WriteLine("overall");
            StringBuilder sb = new StringBuilder();
            overall_minimal_isoforms.ForEach(x => sb.Append(string.Join("/", x.Key) + ','));
            sb = sb.Remove(sb.Length - 1, 1);
            Console.WriteLine(sb.ToString());

            foreach (string subset in subsets)
            {
                List<KeyValuePair<List<string>, List<Peptide>>> subset_minimal_isoforms = null;

                if (ANALYSIS_TYPE == AnalysisType.Independent)
                {
                    // do an independent minimal isoform group analysis for each subset
                    List<KeyValuePair<string, List<Peptide>>> subset_isoforms = new List<KeyValuePair<string, List<Peptide>>>(isoforms.Select(x => new KeyValuePair<string, List<Peptide>>(x.Key, new List<Peptide>(x.Value))));
                    string[] subset_components = subset.Split('_');
                    foreach (Peptide peptide in all_peptides)
                    {
                        bool in_subset = false;
                        foreach (string experiment in peptide.Experiments)
                            if (Array.TrueForAll(subset_components, x => experiment.Contains(x)))
                            {
                                in_subset = true;
                                break;
                            }
                        if (!in_subset)
                            foreach (KeyValuePair<List<string>, List<Peptide>> kvp in subset_minimal_isoforms)
                                foreach (Peptide peptide2 in kvp.Value)
                                    if (peptide2 == peptide)
                                    {
                                        // this peptide isn't found in the current subset; remove it from all isoforms
                                        kvp.Value.Remove(peptide2);
                                        break;
                                    }
                    }
                    subset_isoforms.RemoveAll(x => x.Value.Count == 0);  // remove all isoforms that no longer contain any peptides from this subset
                    subset_minimal_isoforms = FindMinimalProteinGroups(subset_isoforms);
                }
                else if (ANALYSIS_TYPE == AnalysisType.ConserveGroups)
                {
                    // start with the overall set of minimal protein groups and remove peptides not found in this subset; then remove all protein groups that have no remaining peptides or no remaining unique peptides
                    subset_minimal_isoforms = new List<KeyValuePair<List<string>, List<Peptide>>>(overall_minimal_isoforms.Select(x => new KeyValuePair<List<string>, List<Peptide>>(x.Key, new List<Peptide>(x.Value))));
                    string[] subset_components = subset.Split('_');
                    foreach (Peptide peptide in all_peptides)
                    {
                        bool in_subset = false;
                        foreach (string experiment in peptide.Experiments)
                            if (Array.TrueForAll(subset_components, x => experiment.Contains(x)))
                            {
                                in_subset = true;
                                break;
                            }
                        if (!in_subset)
                            foreach (KeyValuePair<List<string>, List<Peptide>> kvp in subset_minimal_isoforms)
                                foreach (Peptide peptide2 in kvp.Value)
                                    if (peptide2 == peptide)
                                    {
                                        // this peptide isn't found in the current subset; remove it from all isoforms
                                        kvp.Value.Remove(peptide2);
                                        break;
                                    }
                    }
                    subset_minimal_isoforms.RemoveAll(x => x.Value.Count == 0);  // remove all isoforms that no longer contain any peptides from this subset
                    if (REMOVE_ISOFORMS_WITHOUT_UNIQUE_PEPTIDES)
                    {
                        // remove all isoforms without unique peptides
                        int i = 0;
                        while (i < subset_minimal_isoforms.Count)
                        {
                            KeyValuePair<List<string>, List<Peptide>> kvp = subset_minimal_isoforms[i];
                            bool has_unique_peptide = false;
                            foreach (Peptide peptide in kvp.Value)
                            {
                                bool peptide_is_unique = true;
                                foreach (KeyValuePair<List<string>, List<Peptide>> kvp2 in subset_minimal_isoforms)
                                {
                                    if (kvp2.Equals(kvp))
                                        continue;
                                    if (kvp2.Value.Contains(peptide))
                                    {
                                        peptide_is_unique = false;
                                        break;
                                    }
                                }
                                if (peptide_is_unique)
                                {
                                    has_unique_peptide = true;
                                    break;
                                }
                            }
                            if (!has_unique_peptide)
                                subset_minimal_isoforms.RemoveAt(i);
                            else
                                i++;
                        }
                    }
                }

                Console.Error.WriteLine(subset + '\t' + subset_minimal_isoforms.Count.ToString());

                Console.WriteLine(subset);
                sb = new StringBuilder();
                subset_minimal_isoforms.ForEach(x => sb.Append(string.Join("/", x.Key) + ','));
                sb = sb.Remove(sb.Length - 1, 1);
                Console.WriteLine(sb.ToString());
            }
        }

        static List<KeyValuePair<List<string>, List<Peptide>>> FindMinimalProteinGroups(ICollection<KeyValuePair<string, List<Peptide>>> maximalProteins)
        {
            List<KeyValuePair<List<string>, List<Peptide>>> minimal_protein_groups = new List<KeyValuePair<List<string>, List<Peptide>>>(maximalProteins.Select(x => new KeyValuePair<List<string>, List<Peptide>>(new List<string>(new[] { x.Key }), new List<Peptide>(x.Value))));

            // merge indistinguishable proteins into groups and remove subset protein groups
            int i = 0;
            while (i < minimal_protein_groups.Count - 1)
            {
                HashSet<int> first_protein_group_peptide_ids = new HashSet<int>(minimal_protein_groups[i].Value.Select(x => x.Id));
                bool remove_first = false;
                int j = i + 1;
                while (j < minimal_protein_groups.Count)
                {
                    HashSet<int> second_protein_group_peptide_ids = new HashSet<int>(minimal_protein_groups[j].Value.Select(x => x.Id));
                    if (second_protein_group_peptide_ids.SetEquals(first_protein_group_peptide_ids))
                    {
                        minimal_protein_groups[i].Key.AddRange(minimal_protein_groups[j].Key);
                        minimal_protein_groups.RemoveAt(j);
                    }
                    else if (second_protein_group_peptide_ids.IsSubsetOf(first_protein_group_peptide_ids))
                        minimal_protein_groups.RemoveAt(j);
                    else if (first_protein_group_peptide_ids.IsSubsetOf(second_protein_group_peptide_ids))
                    {
                        remove_first = true;
                        break;
                    }
                    else
                        j++;
                }
                if (remove_first)
                    minimal_protein_groups.RemoveAt(i);
                else
                    i++;
            }

            // remove subsumable protein groups
            minimal_protein_groups.Sort((l, r) => l.Value.Sum(x => x.Score).CompareTo(r.Value.Sum(x => x.Score)));
            i = 0;
            while (i < minimal_protein_groups.Count)
            {
                HashSet<int> lower_protein_group_peptide_ids = new HashSet<int>(minimal_protein_groups[i].Value.Select(x => x.Id));
                bool remove = false;
                for (int j = 0; j < minimal_protein_groups.Count; j++)
                {
                    if (j != i)
                    {
                        HashSet<int> higher_protein_group_peptide_ids = new HashSet<int>(minimal_protein_groups[j].Value.Select(x => x.Id));
                        lower_protein_group_peptide_ids.ExceptWith(higher_protein_group_peptide_ids);
                        if (lower_protein_group_peptide_ids.Count == 0)
                        {
                            remove = true;
                            break;
                        }
                    }
                }
                if (remove)
                    minimal_protein_groups.RemoveAt(i);
                else
                    i++;
            }

            if (REMOVE_EACH_GROUP_SEQUENTIALLY_TEST)
            {
                HashSet<int> peptide_ids_from_minimal_protein_groups = new HashSet<int>(minimal_protein_groups.SelectMany(x => x.Value.Select(y => y.Id)));
                foreach (KeyValuePair<List<string>, List<Peptide>> protein_group in minimal_protein_groups)
                    if (new HashSet<int>(minimal_protein_groups.Except(new[] { protein_group }).SelectMany(x => x.Value.Select(y => y.Id))).SetEquals(peptide_ids_from_minimal_protein_groups))
                        throw new Exception();
            }

            return minimal_protein_groups;
        }
    }

    class Peptide
    {
        public int Id { get; private set; }
        public double Score { get; private set; }
        public HashSet<string> Experiments { get; private set; }

        public Peptide(int id, double score, HashSet<string> experiments)
        {
            Id = id;
            Score = score;
            Experiments = experiments;
        }
    }
}
