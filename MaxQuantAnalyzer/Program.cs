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
            bool exclude_special_proteins = true;
            if (args.Length > 3)
                if (!bool.TryParse(args[3], out exclude_special_proteins))
                    exclude_special_proteins = true;

            // read in peptides.txt and make a list of all the peptides containg the sequence and the number of times it's found in each experiment, keyed by peptide ID
            Dictionary<int, Tuple<string, Dictionary<string, int>>> peptides = new Dictionary<int, Tuple<string, Dictionary<string, int>>>();
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

                    Dictionary<string, int> psms_per_exp = new Dictionary<string, int>();
                    foreach (KeyValuePair<int, string> kvp in experiment_indexes)
                    {
                        int psms;
                        int.TryParse(fields[kvp.Key], out psms);
                        psms_per_exp.Add(kvp.Value, psms);
                    }

                    peptides.Add(int.Parse(fields[id_index]), new Tuple<string, Dictionary<string, int>>(fields[seq_index], psms_per_exp));
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
                // sequence coverage output
                using (StreamWriter sequence_coverage_out = new StreamWriter(Path.Combine(Path.GetDirectoryName(proteins_filename), Path.GetFileNameWithoutExtension(proteins_filename) + ".sequence_coverage" + Path.GetExtension(proteins_filename))))
                {
                    sequence_coverage_out.AutoFlush = true;

                    // PSMs and unique peptide output
                    using (StreamWriter PSMs_and_unique_peptides_out = new StreamWriter(Path.Combine(Path.GetDirectoryName(proteins_filename), Path.GetFileNameWithoutExtension(proteins_filename) + ".PSMs_and_unique_peptides" + Path.GetExtension(proteins_filename))))
                    {
                        PSMs_and_unique_peptides_out.AutoFlush = true;

                        // non-redundant amino acids output
                        using (StreamWriter nonredundant_amino_acids_out = new StreamWriter(Path.Combine(Path.GetDirectoryName(proteins_filename), Path.GetFileNameWithoutExtension(proteins_filename) + ".non-redundant_amino_acids" + Path.GetExtension(proteins_filename))))
                        {
                            nonredundant_amino_acids_out.AutoFlush = true;

                            string header = proteins_file.ReadLine();
                            string[] header_fields = header.Split('\t');
                            int id_index = Array.IndexOf(header_fields, "id");
                            int length_index = Array.IndexOf(header_fields, "Sequence length");
                            int peptide_ids_index = Array.IndexOf(header_fields, "Peptide IDs");
                            int maj_prot_ids_index = Array.IndexOf(header_fields, "Majority protein IDs");
                            int seq_cov_index = Array.IndexOf(header_fields, "Sequence coverage [%]");
                            int evidence_ids_index = Array.IndexOf(header_fields, "Evidence IDs");
                            int site_only_index = Array.IndexOf(header_fields, "Only identified by site");
                            int reverse_index = Array.IndexOf(header_fields, "Reverse");
                            int contaminant_index = Array.IndexOf(header_fields, "Potential contaminant");
                            Dictionary<string, int> sc_indexes = new Dictionary<string, int>();
                            StringBuilder sequence_coverage_sb = new StringBuilder("id\tProtein ID\t");
                            foreach (string subset in subsets)
                            {
                                int index = Array.IndexOf(header_fields, "Sequence coverage " + subset + " [%]");
                                // if MaxQuant reported sequence coverage for this subset, add a column header to the output and save the column index
                                if (index >= 0)
                                {
                                    sequence_coverage_sb.Append("MaxQuant " + header_fields[index] + '\t');
                                    sc_indexes.Add(subset, index);
                                }
                                sequence_coverage_sb.Append("Sequence coverage " + subset + " [%]\t");
                            }
                            sequence_coverage_sb.Append("MaxQuant Sequence coverage [%]\tSequence coverage [%]");
                            sequence_coverage_out.WriteLine(sequence_coverage_sb.ToString());
                            Dictionary<string, int> distinct_peptides_indexes = new Dictionary<string, int>();
                            StringBuilder PSMs_and_unique_peptides_sb = new StringBuilder("id\tProtein ID\t");
                            foreach (string subset in subsets)
                            {
                                // MaxQuant never reports a number of PSMs on a per-experiment basis so only output our results
                                PSMs_and_unique_peptides_sb.Append("PSMs " + subset + '\t');
                                int index = Array.IndexOf(header_fields, "Peptides " + subset);
                                // if MaxQuant reported a number of distinct peptides for this subset, add a column header to the output and save the column index
                                if (index >= 0)
                                {
                                    PSMs_and_unique_peptides_sb.Append("MaxQuant unique peptides " + subset + '\t');
                                    distinct_peptides_indexes.Add(subset, index);
                                }
                                PSMs_and_unique_peptides_sb.Append("Unique peptides " + subset + '\t');
                            }
                            PSMs_and_unique_peptides_sb.Append("MaxQuant PSMs\tPSMs\tMaxQuant unique peptides\tCDW unique peptides");
                            PSMs_and_unique_peptides_out.WriteLine(PSMs_and_unique_peptides_sb.ToString());
                            nonredundant_amino_acids_out.WriteLine("id\tProtein ID\tTrypsin non-redundant amino acids\tNon-trypsin non-redundant amino acids\tTrypsin and non-trypsin non-redundant amino acids\tTotal non-redundant amino acids");

                            while (!proteins_file.EndOfStream)
                            {
                                string line = proteins_file.ReadLine();
                                string[] fields = line.Split('\t');

                                int id = int.Parse(fields[id_index]);
                                int length = int.Parse(fields[length_index]);
                                int[] peptide_ids = Array.ConvertAll(fields[peptide_ids_index].Split(';'), x => int.Parse(x));
                                string[] maj_prot_ids = fields[maj_prot_ids_index].Split(';');
                                string maj_prot_id = maj_prot_ids[0];

                                if (exclude_special_proteins && (!string.IsNullOrWhiteSpace(fields[site_only_index]) || !string.IsNullOrWhiteSpace(fields[reverse_index]) || !string.IsNullOrWhiteSpace(fields[contaminant_index])))
                                    continue;

                                // get the protein sequence of the first major protein (ignore if reversed/decoy or contaminant)
                                string maj_prot_seq = null;
                                if (!maj_prot_id.StartsWith("REV") && !maj_prot_id.StartsWith("CON"))
                                    maj_prot_seq = protein_sequences[maj_prot_id];

                                // for each subset of experiments, calculate sequence coverage and number of PSMs and distinct peptides for each protein group
                                sequence_coverage_sb = new StringBuilder(fields[id_index] + '\t' + maj_prot_id + '\t');
                                PSMs_and_unique_peptides_sb = new StringBuilder(fields[id_index] + '\t' + maj_prot_id + '\t');
                                foreach (string subset in subsets)
                                {
                                    string[] subset_components = subset.Split('_');

                                    double seq_cov = double.NaN;
                                    int psms = 0;
                                    int distinct_peptides = 0;
                                    if (maj_prot_seq != null)
                                    {
                                        // for every peptide associated with the protein, figure out if it was found in the current subset of runs
                                        // if so, find every instance of it in the parent protein sequence
                                        // mark those residues as covered
                                        // after this is done for all peptides, calculate sequence coverage
                                        HashSet<int> residues = new HashSet<int>();
                                        foreach (int peptide_id in peptide_ids)
                                        {
                                            Tuple<string, Dictionary<string, int>> peptide;
                                            if (!peptides.TryGetValue(peptide_id, out peptide))
                                                continue;

                                            bool in_subset = false;
                                            foreach (KeyValuePair<string, int> kvp in peptide.Item2)
                                                if (kvp.Value > 0 && Array.TrueForAll(subset_components, x => kvp.Key.Contains(x)))
                                                {
                                                    in_subset = true;
                                                    psms += kvp.Value;
                                                }

                                            if (in_subset)
                                            {
                                                MatchCollection matches = Regex.Matches(maj_prot_seq, peptide.Item1);
                                                if (matches.Count == 0)
                                                    throw new Exception("Peptide not found in parent protein.");

                                                foreach (Match match in matches)
                                                    for (int r = match.Index + 1; r < match.Index + 1 + match.Length; r++)
                                                        residues.Add(r);

                                                distinct_peptides++;
                                            }
                                        }
                                        seq_cov = (double)residues.Count / length * 100;
                                    }
                                    int index;
                                    if (sc_indexes.TryGetValue(subset, out index))
                                        sequence_coverage_sb.Append(fields[index] + '\t');
                                    sequence_coverage_sb.Append(seq_cov.ToString("F1") + '\t');
                                    PSMs_and_unique_peptides_sb.Append(psms.ToString() + '\t');
                                    if (distinct_peptides_indexes.TryGetValue(subset, out index))
                                        PSMs_and_unique_peptides_sb.Append(fields[index] + '\t');
                                    PSMs_and_unique_peptides_sb.Append(distinct_peptides.ToString() + '\t');
                                }

                                // calculate overall sequence coverage, number of PSMs, number of distinct peptides, number of tryptic/non-tryptic/intersection(tryptic, non-tryptic) non-redundant amino acids for each protein group
                                double overall_seq_cov = double.NaN;
                                int overall_psms = 0;
                                int overall_distinct_peptides = 0;
                                HashSet<int> total_residues = new HashSet<int>();
                                HashSet<int> tryptic_residues = new HashSet<int>();
                                HashSet<int> non_tryptic_residues = new HashSet<int>();
                                if (maj_prot_seq != null)
                                {
                                    foreach (int peptide_id in peptide_ids)
                                    {
                                        Tuple<string, Dictionary<string, int>> peptide;
                                        if (!peptides.TryGetValue(peptide_id, out peptide))
                                            continue;

                                        bool found_with_trypsin = false;
                                        bool found_without_trypsin = false;
                                        foreach (KeyValuePair<string, int> kvp in peptide.Item2)
                                        {
                                            if (kvp.Value > 0)
                                            {
                                                if (kvp.Key.Contains("Trypsin"))
                                                    found_with_trypsin = true;
                                                else
                                                    found_without_trypsin = true;
                                            }
                                            if (found_with_trypsin && found_without_trypsin)
                                                break;
                                        }

                                        MatchCollection matches = Regex.Matches(maj_prot_seq, peptide.Item1);
                                        if (matches.Count == 0)
                                            throw new Exception("Peptide not found in parent protein.");

                                        foreach (Match match in matches)
                                            for (int r = match.Index + 1; r < match.Index + 1 + match.Length; r++)
                                            {
                                                total_residues.Add(r);
                                                if (found_with_trypsin)
                                                    tryptic_residues.Add(r);
                                                if (found_without_trypsin)
                                                    non_tryptic_residues.Add(r);
                                            }

                                        foreach (KeyValuePair<string, int> kvp in peptide.Item2)
                                            overall_psms += kvp.Value;
                                    }
                                    overall_seq_cov = (double)total_residues.Count / length * 100;
                                    overall_distinct_peptides += peptide_ids.Length;
                                }
                                sequence_coverage_sb.Append(fields[seq_cov_index] + '\t' + overall_seq_cov.ToString("F1"));
                                PSMs_and_unique_peptides_sb.Append(fields[evidence_ids_index].Split(';').Length.ToString() + '\t' + overall_psms.ToString() + '\t' + peptide_ids.Length.ToString() + '\t' + overall_distinct_peptides.ToString());
                                sequence_coverage_out.WriteLine(sequence_coverage_sb.ToString());
                                PSMs_and_unique_peptides_out.WriteLine(PSMs_and_unique_peptides_sb.ToString());
                                HashSet<int> tryptic_and_non_tryptic_residues = new HashSet<int>(tryptic_residues);
                                tryptic_and_non_tryptic_residues.IntersectWith(non_tryptic_residues);
                                nonredundant_amino_acids_out.WriteLine(fields[id_index] + '\t' + maj_prot_id + '\t' + tryptic_residues.Count.ToString() + '\t' + non_tryptic_residues.Count.ToString() + '\t' + tryptic_and_non_tryptic_residues.Count.ToString() + '\t' + total_residues.Count.ToString());
                            }
                        }
                    }
                }
            }
        }
    }
}
