using DotNetWikiBot;
using System;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace wikisharp
{
    class Program
    {
        static readonly Regex LINKED_PIPED_SECTION_HEADING = new Regex(@"(={1,6})\[\[(.+)\|(.+?)\]\](={1,6})");
        static readonly Regex LINKED_UNPIPED_SECTION_HEADING = new Regex(@"(={1,6})\[\[(.+?)\]\](={1,6})");

        static void Main(string[] args)
        {
            Console.Write("Username: ");
            string username = Console.ReadLine();
            Console.Write("Password: ");
            // from https://stackoverflow.com/a/3404522/60067
            string password = string.Empty;
            ConsoleKey key;
            do
            {
                var keyInfo = Console.ReadKey(intercept: true);
                key = keyInfo.Key;
                if (key == ConsoleKey.Backspace && password.Length > 0)
                {
                    Console.Write("\b \b");
                    password = password.Substring(0, password.Length - 1);
                }
                else if (!char.IsControl(keyInfo.KeyChar))
                {
                    Console.Write("*");
                    password += keyInfo.KeyChar;
                }
            } while (key != ConsoleKey.Enter);

            Site s = new Site("https://en.wikipedia.org", username, password);
            PageList pl = new PageList(s);
            pl.FillFromCategory("American Civil War orders of battle");

            foreach (Page p in pl)
            {
                p.Load();
                string wikitext = p.text;
                string modified_wt = Modify(wikitext);
                if (modified_wt != wikitext)
                {
                    p.Save(modified_wt, "[[MOS:NOSECTIONLINKS]]", false);
                    Console.WriteLine(p.title + " modified.");
                }
            }
        }

        static string Modify(string original)
        {
            StringBuilder modified_sb = new StringBuilder(original.Length);
            using (StringReader sr = new StringReader(original))
            {
                while (sr.Peek() != -1)
                {
                    string line = sr.ReadLine();
                    Match piped_match = LINKED_PIPED_SECTION_HEADING.Match(line);
                    if (piped_match.Success)
                    {
                        modified_sb.AppendLine(piped_match.Groups[1].Value + piped_match.Groups[3].Value + piped_match.Groups[4].Value);
                        modified_sb.AppendLine("{{Main|" + piped_match.Groups[2].Value + "}}");
                    }
                    else
                    {
                        Match unpiped_match = LINKED_UNPIPED_SECTION_HEADING.Match(line);
                        if (unpiped_match.Success)
                        {
                            modified_sb.AppendLine(unpiped_match.Groups[1].Value + unpiped_match.Groups[2].Value + unpiped_match.Groups[3].Value);
                            modified_sb.AppendLine("{{Main|" + unpiped_match.Groups[2].Value + "}}");
                        }
                        else
                            modified_sb.AppendLine(line);
                    }
                }
            }
            return modified_sb.ToString().Replace("\r\n", "\n").Trim();
        }
    }
}
