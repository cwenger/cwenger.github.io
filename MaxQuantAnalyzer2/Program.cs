using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MaxQuantAnalyzer2
{
    class Program
    {
        static void Main(string[] args)
        {
            string peptides_filename = args[0];
            bool exclude_nomsms_peptides;
            if (args.Length >= 2)
                bool.TryParse(args[1], out exclude_nomsms_peptides);
            else
                exclude_nomsms_peptides = false;

            HashSet<int> all_peptide_ids = new HashSet<int>();
            Dictionary<string, List<Peptide>> isoforms_dict = new Dictionary<string, List<Peptide>>();
            using (StreamReader peptides_file = new StreamReader(peptides_filename))
            {
                string header = peptides_file.ReadLine();
                string[] header_fields = header.Split('\t');
                int id_index = Array.IndexOf(header_fields, "id");
                int best_msms_index = Array.IndexOf(header_fields, "Best MS/MS");
                int score_index = Array.IndexOf(header_fields, "Score");
                int isoforms_index = Array.IndexOf(header_fields, "isoforms(+)");

                while (!peptides_file.EndOfStream)
                {
                    string line = peptides_file.ReadLine();
                    string[] fields = line.Split('\t');

                    if (exclude_nomsms_peptides && string.IsNullOrWhiteSpace(fields[best_msms_index]))
                        continue;

                    string isoforms_field = fields[isoforms_index];
                    if (string.IsNullOrWhiteSpace(isoforms_field))
                        continue;
                    int id = int.Parse(fields[id_index]);
                    all_peptide_ids.Add(id);
                    Peptide peptide = new Peptide(id, double.Parse(fields[score_index]));
                    foreach (string isoform in isoforms_field.Split(';'))
                    {
                        List<Peptide> peptides;
                        if (!isoforms_dict.TryGetValue(isoform, out peptides))
                        {
                            peptides = new List<Peptide>();
                            peptides.Add(peptide);
                            isoforms_dict.Add(isoform, peptides);
                        }
                        else
                            peptides.Add(peptide);
                    }
                }
            }

            List<KeyValuePair<List<string>, List<Peptide>>> isoforms = new List<KeyValuePair<List<string>, List<Peptide>>>(isoforms_dict.Select(x => new KeyValuePair<List<string>, List<Peptide>>(new List<string>(new[] { x.Key }), x.Value)));

            // merge indistinguishable isoforms and remove subset isoforms
            //isoforms.Sort((l, r) => -(l.Value.Sum(x => x.Score).CompareTo(r.Value.Sum(x => x.Score))));
            int i = 0;
            while (i < isoforms.Count)
            {
                HashSet<int> first_isoform_peptide_ids = new HashSet<int>(isoforms[i].Value.Select(x => x.Id));
                bool remove_first = false;
                int j = 0;
                while (j < isoforms.Count)
                {
                    if (j != i)
                    {
                        HashSet<int> second_isoform_peptide_ids = new HashSet<int>(isoforms[j].Value.Select(x => x.Id));
                        if (second_isoform_peptide_ids.SetEquals(first_isoform_peptide_ids))
                        {
                            isoforms[i].Key.AddRange(isoforms[j].Key);
                            isoforms.RemoveAt(j);
                        }
                        else if (second_isoform_peptide_ids.IsSubsetOf(first_isoform_peptide_ids))
                            isoforms.RemoveAt(j);
                        else if (first_isoform_peptide_ids.IsSubsetOf(second_isoform_peptide_ids))
                        {
                            remove_first = true;
                            break;
                        }
                        else
                            j++;
                    }
                    else
                        j++;
                }
                if (remove_first)
                    isoforms.RemoveAt(i);
                else
                    i++;
            }

            // remove subsumable isoforms
            isoforms.Sort((l, r) => l.Value.Sum(x => x.Score).CompareTo(r.Value.Sum(x => x.Score)));
            i = 0;
            while (i < isoforms.Count)
            {
                HashSet<int> lower_isoform_peptide_ids = new HashSet<int>(isoforms[i].Value.Select(x => x.Id));
                bool remove = false;
                for (int j = 0; j < isoforms.Count; j++)
                {
                    if (j != i)
                    {
                        HashSet<int> higher_isoform_peptide_ids = new HashSet<int>(isoforms[j].Value.Select(x => x.Id));
                        lower_isoform_peptide_ids.ExceptWith(higher_isoform_peptide_ids);
                        if (lower_isoform_peptide_ids.Count == 0)
                        {
                            remove = true;
                            break;
                        }
                    }
                }
                if (remove)
                    isoforms.RemoveAt(i);
                else
                    i++;
            }

            HashSet<int> peptide_ids_from_parsimony_isoforms = new HashSet<int>(isoforms.SelectMany(x => x.Value.Select(y => y.Id)));
            if (!all_peptide_ids.SetEquals(peptide_ids_from_parsimony_isoforms))
                throw new Exception();

            //foreach (KeyValuePair<string, List<Peptide>> isoform in isoforms)
            //    if (new HashSet<int>(isoforms.Except(new[] { isoform }).SelectMany(x => x.Value.Select(y => y.Id))).SetEquals(peptide_ids_from_parsimony_isoforms))
            //        throw new Exception();
        }
    }

    public class Peptide
    {
        public int Id { get; private set; }
        public double Score { get; private set; }

        public Peptide(int id, double score)
        {
            Id = id;
            Score = score;
        }
    }
}
