package model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TmpFileManager 
{
	private TreeSet<Read> reads = new TreeSet<Read>();
	private HashMap<String, Read> topReads = new HashMap<String, Read>();
	private HashMap<String, BufferedReader> whichReader = new HashMap<String, BufferedReader>();
	private ArrayList<SamReader> samReaders = new ArrayList<SamReader>();
	private HashMap<String, SAMRecordIterator> whichReaderSAM = new HashMap<String, SAMRecordIterator>();
	
	public TmpFileManager(String type, int nbTmp) throws Exception
	{
		switch(type)
		{
			case "bam":
			case "fastq":
				for(int i = 1; i <= nbTmp; i++)
				{
					BufferedReader br = new BufferedReader(new FileReader(Parameters.tmpFolder + type + "." + i + ".tmp"));
					Read r = parse(type, br);
					reads.add(r);
					topReads.put(r.name, r);
					whichReader.put(r.name, br);
				}
				break;
			case "sam":
				for(int i = 1; i <= nbTmp; i++)
				{
					SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
					SamReader samReader = samReaderFactory.open(new File(Parameters.tmpFolder + type + "." + i + ".tmp"));
					SAMRecordIterator it = samReader.iterator();
					Read r = parse(it);
					reads.add(r);
					topReads.put(r.name, r);
					whichReaderSAM.put(r.name, it);
					samReaders.add(samReader);
				}
		}
	}
	
	public static void createTmpFile(String filename, TreeSet<String> lines) throws IOException
	{
		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
		for(String l:lines) bw.write(l + "\n");
		bw.close();
	}
	
	public static void createTmpFileB(SAMFileHeader header, String filename, TreeSet<ComparableBAMRecord> lines) throws IOException
	{
		SAMFileWriterFactory samWriterFactory = new SAMFileWriterFactory();
		SAMFileWriter samWriter = samWriterFactory.makeBAMWriter(header, true, new File(filename));
		for(ComparableBAMRecord l:lines) samWriter.addAlignment(l.record);
		samWriter.close();
	}
		
	private Read parse(String type, BufferedReader br) throws IOException
	{
		String line = br.readLine();
		if(line == null) return null;
		String[] tokens = line.split("\t");
		Read r = new Read();
		r.name = tokens[0];
		switch(type)
		{
		case "bam":
			r.gene = tokens[1];
			break;
		case "fastq":
			r.barcode = tokens[1];
			r.UMI = tokens[2];
			if(tokens.length == 5)
			{
				r.qualityBC = tokens[3];
				r.qualityUMI = tokens[4];
			}
		}
		return r;
	}
	
	private Read parse(SAMRecordIterator it) throws IOException
	{
		if(!it.hasNext()) return null;
		SAMRecord record = it.next();
		Read r = new Read();
		r.name = record.getReadName();
		r.samRecord = record;
		r.gene = (String)r.samRecord.getAttribute("BB"); // Gene info
		return r;
	}
	
	public static void delete(String type, int nbTmp) throws Exception
	{
		for(int i = 1; i <= nbTmp; i++) new File(Parameters.tmpFolder + type + "." + i + ".tmp").delete();
	}
	
	public Read getFirstRead() throws Exception // Always called for fastq file
	{
		Read r = reads.pollFirst();
		if(r == null) return null;
		// Get next read from the corresponding reader
		BufferedReader br = whichReader.get(r.name);
		Read r2 = parse("fastq", br);
		if(r2 != null) // If null the readName is removed and nothing is added back
		{
			reads.add(r2);
			topReads.put(r2.name, r2);
			whichReader.put(r2.name, br);
		}
		else br.close(); // I close the reader when it's over
		// Remove from top
		topReads.remove(r.name);
		whichReader.remove(r.name);
		return r;
	}
	
	public Read getReadB(String readName) throws Exception
	{
		Read r = topReads.get(readName);
		if(r != null) // This read is indeed the top of one tmp file
		{
			BufferedReader br = whichReader.get(readName);
			Read r2 = parse("bam", br);
			if(r2 != null) // If null the readName is removed and nothing is added back
			{
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
	
	public Read getReadS(String readName) throws Exception
	{
		Read r = topReads.get(readName);
		if(r != null) // This read is indeed the top of one tmp file
		{
			SAMRecordIterator it = whichReaderSAM.get(readName);
			Read r2 = parse(it);
			if(r2 != null) // If null the readName is removed and nothing is added back
			{
				reads.add(r2);
				topReads.put(r2.name, r2);
				whichReaderSAM.put(r2.name, it);
			}
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
	
	public void closeSAMReaders() throws IOException
	{
		for(SamReader reader:samReaders) reader.close();
	}
}
