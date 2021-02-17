package model;

import java.io.BufferedReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;
import tools.Utils;

public class GTF 
{
	private static HashMap<String, IntervalTree> forest = null;
	
	public static void readGTF() throws Exception
	{
		System.out.println("\nReading GTF file provided: " + Parameters.inputGTFFile.getAbsolutePath());
		Parameters.geneIndex = new HashMap<String, Integer>(); // Fill this as we go
		Parameters.mappingGeneIdGeneName = new HashMap<String, HashSet<String>>(); // For filling final matrices
		BufferedReader br = Utils.readGTF(Parameters.inputGTFFile);
		String line = br.readLine();
		int line_number = 0;
		forest = new HashMap<>();
		int nbExons = 0;
		int nbGenes = 0;
		while(line != null)
		{
			line_number++;
	    	if(!line.startsWith("#") && !line.trim().equals(""))
    		{
				String[] tokens = line.split("\t");
				// Parse Line
				String[] params = tokens[8].split(";");
				long start = Long.parseLong(tokens[3]);
				long end = Long.parseLong(tokens[4]);
				String gene_name = null;
				String gene_id = null;
				String chr = tokens[0];
				String type = tokens[2];
				boolean strand = tokens[6].equals("+");
				for(String param:params) 
				{
					String[] values = param.trim().split("\\s+");
					if(values.length >= 2)
					{
						values[1] = values[1].replaceAll("\"", "");
						if(values[0].equals("gene_name")) gene_name = values[1];
						else if(values[0].equals("gene_id")) gene_id = values[1];
					}
				}
				if(gene_name == null) gene_name = gene_id;
				if(gene_id == null) System.err.println("[WARNING] l." + line_number + "\t: No gene_id. This entry is ignored.");
				else
				{
					HashSet<String> names = Parameters.mappingGeneIdGeneName.get(gene_id);
					if(names == null) // First time we see this gene_id
					{
						names = new HashSet<String>();
						Parameters.geneIndex.put(gene_id, nbGenes);
						nbGenes++;
					}
					names.add(gene_name);
					Parameters.mappingGeneIdGeneName.put(gene_id, names);
					// Which type is it?
					if(type.equals("exon")) 
					{
						nbExons++;
						IntervalTree tree = forest.get(chr);
						if(tree == null) tree = new IntervalTree();
						tree.insert(new IntervalLabelled((int)start, (int)end, gene_id, strand));
						forest.put(chr, tree);
					}
				}
			}
			line = br.readLine();
		}
		br.close();

		if(nbGenes == 0) {
			System.err.println("We couldn't parse the GTF file. Please report this problem if the GTF is in standard format. Or use another GTF from another source.");
			System.exit(-1);
		}
		
		System.out.println(nbExons + " 'exons' are annotating " + nbGenes + " unique gene_ids in the provided GTF file.");
		
		Parameters.geneIndex.put("__alignment_not_unique", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__no_feature", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__ambiguous", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__too_low_aQual", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__not_aligned", Parameters.geneIndex.size());
	}
	
	public static HashSet<String> findOverlappingGenes(String chr, int start, int end, boolean readNegativeStrandFlag)
	{
		HashSet<String> result = new HashSet<>();
		IntervalTree itree = forest.get(chr);
		if(itree != null)
		{
			List<Interval> found = itree.findOverlapping(new Interval(start, end));
			for(Interval i:found) 
			{
				IntervalLabelled g = (IntervalLabelled)i;
				if(Parameters.stranded == Strand.NO) result.add(g.gene);
				else if(Parameters.stranded == Strand.REVERSE && g.readNegativeStrandFlag == readNegativeStrandFlag) result.add(g.gene);
				else if(Parameters.stranded == Strand.YES && g.readNegativeStrandFlag != readNegativeStrandFlag) result.add(g.gene);
			}
		}
		return result;
	}
}
