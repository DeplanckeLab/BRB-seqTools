package model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.TreeSet;

public class TmpFileManager 
{
	private TreeSet<Read> reads = new TreeSet<Read>();
	private HashMap<String, Read> topReads = new HashMap<String, Read>();
	private HashMap<String, BufferedReader> whichReader = new HashMap<String, BufferedReader>();
	
	public TmpFileManager(String type, int nbTmp) throws Exception
	{
		for(int i = 1; i <= nbTmp; i++)
		{
			BufferedReader br = new BufferedReader(new FileReader(Parameters.tmpFolder + type + "." + i + ".tmp"));
			Read r = parse(br.readLine());
			reads.add(r);
			topReads.put(r.name, r);
			whichReader.put(r.name, br);
		}
	}
	
	private Read parse(String line)
	{
		String[] tokens = line.split("\t");
		Read r = new Read();
		r.name = tokens[0];
		if(tokens.length > 2)
		{
			r.barcode = tokens[1];
			r.UMI = tokens[2];
		}
		else r.gene = tokens[1];
		return r;
	}
	
	public static void delete(String type, int nbTmp) throws Exception
	{
		for(int i = 1; i <= nbTmp; i++) new File(Parameters.tmpFolder + type + "." + i + ".tmp").delete();
	}
	
	public Read getFirstRead() throws Exception
	{
		Read r = reads.pollFirst();
		if(r == null) return null;
		// Get next read from the corresponding reader
		BufferedReader br = whichReader.get(r.name);
		String line = br.readLine();  
		if(line != null) // If null the readName is removed and nothing is added back
		{
			Read r2 = parse(line);
			reads.add(r2);
			topReads.put(r2.name, r2);
			whichReader.put(r2.name, br);
		} else br.close();
		// Remove from top
		topReads.remove(r.name);
		whichReader.remove(r.name);
		return r;
	}
	
	public Read getRead(String readName) throws Exception
	{
		Read r = topReads.get(readName);
		if(r != null) // This read is indeed the top of one tmp file
		{
			BufferedReader br = whichReader.get(readName);
			String line = br.readLine();  // Get next read from this reader
			if(line != null) // If null the readName is removed and nothing is added back
			{
				Read r2 = parse(line);
				reads.add(r2);
				topReads.put(r2.name, r2);
				whichReader.put(r2.name, br);
			}
			else br.close(); // I close the reader when it's over
			// Remove from top
			reads.remove(r);
			topReads.remove(readName);
			whichReader.remove(readName);
			return r;
		}
		else // this is not a read mapped
		{
			Parameters.unmapped++;
			return null; 
		}
	}
}
