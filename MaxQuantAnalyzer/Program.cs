using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace MaxQuantAnalyzer
{
    class Program
    {
        static readonly Regex EXPERIMENT_MATCH = new Regex(@"(.+)_(.+)_(.+)");
        static readonly Regex SP_TR_REGEX = new Regex(@"(?:(?:sp|tr)\|(.+)\|)");
        static readonly Regex OTHER_REGEX = new Regex(@"(.+?) ");

        static void Main(string[] args)
        {
            string proteins_filename = args[0];
            string peptides_filename = args[1];
            string[] protein_sequence_filenames = args[2].Split(';');

            // read in peptides.txt and make a list of all the peptides containg the sequence and the experiments it's found in, keyed by peptide ID
            Dictionary<int, Tuple<string, HashSet<string>>> peptides = new Dictionary<int, Tuple<string, HashSet<string>>>();
            HashSet<string> subsets = new HashSet<string>();
            using (StreamReader peptides_file = new StreamReader(peptides_filename))
            {
                string header = peptides_file.ReadLine();
                string[] header_fields = header.Split('\t');
                int id_index = Array.IndexOf(header_fields, "id");
                int seq_index = Array.IndexOf(header_fields, "Sequence");
                Dictionary<int, string> experiment_indexes = new Dictionary<int, string>();
                for (int h = 0; h < header_fields.Length; h++)
                {
                    // make a list of each possible subset of runs
                    // also create a list of all experiments keyed by their column index
                    if (header_fields[h].StartsWith("Experiment"))
                    {
                        string experiment = header_fields[h].Substring("Experiment ".Length);
                        Match match = EXPERIMENT_MATCH.Match(experiment);
                        subsets.Add(match.Groups[1].Value);
                        subsets.Add(match.Groups[2].Value);
                        subsets.Add(match.Groups[3].Value);
                        subsets.Add(match.Groups[1].Value + '_' + match.Groups[2].Value);
                        subsets.Add(match.Groups[2].Value + '_' + match.Groups[3].Value);
                        subsets.Add(match.Groups[1].Value + '_' + match.Groups[3].Value);
                        subsets.Add(experiment);
                        experiment_indexes.Add(h, experiment);
                    }
                }

                while (!peptides_file.EndOfStream)
                {
                    string line = peptides_file.ReadLine();
                    string[] fields = line.Split('\t');

                    HashSet<string> experiments_found_in = new HashSet<string>();
                    foreach (KeyValuePair<int, string> kvp in experiment_indexes)
                        if (!string.IsNullOrWhiteSpace(fields[kvp.Key]))
                            experiments_found_in.Add(kvp.Value);

                    peptides.Add(int.Parse(fields[id_index]), new Tuple<string, HashSet<string>>(fields[seq_index], experiments_found_in));
                }
            }

            // read in all protein sequences, keyed by identifier
            Dictionary<string, string> protein_sequences = new Dictionary<string, string>();
            foreach (string protein_sequence_filename in protein_sequence_filenames)
            {
                using (StreamReader fasta = new StreamReader(protein_sequence_filename))
                {
                    string description = null;
                    string sequence = null;

                    while (true)
                    {
                        if (fasta.EndOfStream || fasta.Peek() == '>')
                        {
                            if (!string.IsNullOrWhiteSpace(description) && !string.IsNullOrWhiteSpace(sequence))
                            {
                                Match sr_tr_match = SP_TR_REGEX.Match(description);
                                if (sr_tr_match.Success)
                                    description = sr_tr_match.Groups[1].Value;
                                else
                                {
                                    Match other_match = OTHER_REGEX.Match(description);
                                    if (other_match.Success)
                                        description = other_match.Groups[1].Value;
                                }
                                protein_sequences.Add(description, sequence);
                                sequence = null;
                            }
                            if (fasta.EndOfStream)
                                break;
                            string line = fasta.ReadLine();
                            description = line.Substring(1);
                        }
                        else
                        {
                            string line = fasta.ReadLine();
                            sequence += line;
                        }
                    }
                }
            }

            using (StreamReader proteins_file = new StreamReader(proteins_filename))
            {
                using (StreamWriter proteins_out_file = new StreamWriter(Path.Combine(Path.GetDirectoryName(proteins_filename), Path.GetFileNameWithoutExtension(proteins_filename) + ".out" + Path.GetExtension(proteins_filename))))
                {
                    proteins_out_file.AutoFlush = true;

                    string header = proteins_file.ReadLine();
                    string[] header_fields = header.Split('\t');
                    int id_index = Array.IndexOf(header_fields, "id");
                    int length_index = Array.IndexOf(header_fields, "Sequence length");
                    int peptide_ids_index = Array.IndexOf(header_fields, "Peptide IDs");
                    int maj_prot_ids_index = Array.IndexOf(header_fields, "Majority protein IDs");
                    int seq_cov_index = Array.IndexOf(header_fields, "Sequence coverage [%]");
                    Dictionary<string, int> sc_indexes = new Dictionary<string, int>();
                    StringBuilder sb = new StringBuilder("id\tProtein ID\t");
                    foreach (string subset in subsets)
                    {
                        int index = Array.IndexOf(header_fields, "Sequence coverage " + subset + " [%]");
                        if (index >= 0)
                        {
                            sb.Append("MaxQuant " + header_fields[index] + '\t');
                            sc_indexes.Add(subset, index);
                        }
                        sb.Append("CDW Sequence coverage " + subset + " [%]\t");
                    }
                    sb.Append("MaxQuant Sequence coverage [%]\tCDW Sequence coverage [%]");
                    proteins_out_file.WriteLine(sb.ToString());

                    while (!proteins_file.EndOfStream)
                    {
                        string line = proteins_file.ReadLine();
                        string[] fields = line.Split('\t');

                        int id = int.Parse(fields[id_index]);
                        int length = int.Parse(fields[length_index]);
                        int[] peptide_ids = Array.ConvertAll(fields[peptide_ids_index].Split(';'), x => int.Parse(x));
                        string[] maj_prot_ids = fields[maj_prot_ids_index].Split(';');
                        string maj_prot_id = maj_prot_ids[0];

                        // get the protein sequence of the first major protein (ignore if reversed/decoy or contaminant)
                        string maj_prot_seq = null;
                        if (!maj_prot_id.StartsWith("REV") && !maj_prot_id.StartsWith("CON"))
                            maj_prot_seq = protein_sequences[maj_prot_id];

                        sb = new StringBuilder(fields[id_index] + '\t' + maj_prot_id + '\t');
                        foreach (string subset in subsets)
                        {
                            string[] subset_components = subset.Split('_');

                            double seq_cov = double.NaN;
                            if (maj_prot_seq != null)
                            {
                                // for every peptide associated with the protein, figure out if it was found in the current subset of runs
                                // if so, find every instance of it in the parent protein sequence
                                // mark those residues as covered
                                // after this is done for all peptides, calculate sequence coverage
                                HashSet<int> residues = new HashSet<int>();
                                foreach (int peptide_id in peptide_ids)
                                {
                                    Tuple<string, HashSet<string>> peptide = peptides[peptide_id];

                                    bool in_subset = false;
                                    foreach (string experiments_found_in in peptide.Item2)
                                    {
                                        if (Array.TrueForAll(subset_components, x => experiments_found_in.Contains(x)))
                                        {
                                            in_subset = true;
                                            break;
                                        }
                                    }

                                    if (in_subset)
                                    {
                                        MatchCollection matches = Regex.Matches(maj_prot_seq, peptide.Item1);
                                        if (matches.Count == 0)
                                            throw new Exception("Peptide not found in parent protein.");

                                        foreach (Match match in matches)
                                            for (int r = match.Index + 1; r < match.Index + 1 + match.Length; r++)
                                                residues.Add(r);
                                    }
                                }
                                seq_cov = (double)residues.Count / length * 100;
                            }
                            int index;
                            if (sc_indexes.TryGetValue(subset, out index))
                                sb.Append(fields[index] + '\t' + seq_cov.ToString("F1") + '\t');
                            else
                                sb.Append(seq_cov.ToString("F1") + '\t');
                        }

                        // calculate overall sequence
                        double overall_seq_cov = double.NaN;
                        if (maj_prot_seq != null)
                        {
                            HashSet<int> residues = new HashSet<int>();
                            foreach (int peptide_id in peptide_ids)
                            {
                                Tuple<string, HashSet<string>> peptide = peptides[peptide_id];

                                MatchCollection matches = Regex.Matches(maj_prot_seq, peptide.Item1);
                                if (matches.Count == 0)
                                    throw new Exception();

                                foreach (Match match in matches)
                                    for (int r = match.Index + 1; r < match.Index + 1 + match.Length; r++)
                                        residues.Add(r);
                            }
                            overall_seq_cov = (double)residues.Count / length * 100;
                        }
                        sb.Append(fields[seq_cov_index] + '\t' + overall_seq_cov.ToString("F1"));
                        proteins_out_file.WriteLine(sb.ToString());
                    }
                }
            }
        }
    }
}
