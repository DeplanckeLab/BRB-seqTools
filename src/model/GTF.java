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
		Parameters.mappingGeneIdGeneName = new HashMap<String, String>(); // For filling final matrices
		BufferedReader br = Utils.readGTF(Parameters.inputGTFFile);
		String line = br.readLine();
		forest = new HashMap<>();
		int nbExons = 0;
		int nbGenes = 0;
		HashSet<String> uniqueGeneId = new HashSet<String>();
		while(line != null)
		{
	    	if(!line.startsWith("#"))
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
				for(String param:params) 
				{
					String value = param.substring(param.indexOf("\"")+1, param.lastIndexOf("\""));
					if(param.contains("gene_name")) gene_name = value;
					if(param.contains("gene_id")) gene_id = value;
				}
				if(gene_name == null) gene_name = gene_id;
				// Which type is it?
				if(type.equals("exon")) 
				{
					nbExons++;
					IntervalTree tree = forest.get(chr);
					if(tree == null) tree = new IntervalTree();
					uniqueGeneId.add(gene_id);
					tree.insert(new IntervalLabelled((int)start, (int)end, gene_id));
					forest.put(chr, tree);
				}
				else if(type.equals("gene"))
				{
					Parameters.geneIndex.put(gene_id, nbGenes);
					Parameters.mappingGeneIdGeneName.put(gene_id, gene_name);
					nbGenes++;
				}
			}
			line = br.readLine();
		}
		br.close();
		
		Parameters.geneIndex.put("__alignment_not_unique", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__no_feature", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__ambiguous", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__too_low_aQual", Parameters.geneIndex.size());
		Parameters.geneIndex.put("__not_aligned", Parameters.geneIndex.size());
		
		System.out.println(nbExons + " 'exons' are annotating " + uniqueGeneId.size() + " unique genes in the provided GTF file. In total " + nbGenes + " 'gene' annotations are found in the GTF file.");
	}
	
	public static HashSet<String> findOverlappingGenes(String chr, int start, int end)
	{
		HashSet<String> result = new HashSet<>();
		IntervalTree itree = forest.get(chr);
		if(itree != null)
		{
			List<Interval> found = itree.findOverlapping(new Interval(start, end));
			for(Interval i:found) result.add(((IntervalLabelled)i).gene);
		}
		return result;
	}
}
