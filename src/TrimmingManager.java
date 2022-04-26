import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import model.Excel;
import model.Parameters;
import tools.Utils;

public class TrimmingManager 
{
	public static final int MAXREADSIZE = 1000;
	public static int[] contaminated = new int[MAXREADSIZE];
	public static int[] polyADist_before = new int[MAXREADSIZE];
	public static int[] polyADist_after = new int[MAXREADSIZE];
	public static String currentFilename = "";

	private static Result searchFor(String search, String in, int oneMismatchAllowedEveryXXbp, int minOverlap)
	{
		Result res = new Result();
		// The "search" object if my BRBseq construct. I will slide from right to left (minimum 10 overlaps) to find it in the "in" String (the difficulty is that it can be cut).
		for(int start = in.length() - minOverlap; start >= -search.length() + minOverlap; start--) 
		{
			int lE = Math.min(start + search.length(), in.length()); // end of in
			int lS = Math.max(0, lE - in.length() + start); // beginning of in
			int allowedMismatch = (lE - lS) / oneMismatchAllowedEveryXXbp; // 1 mismatch allowed every 10 bp for e.g.
			int mismatch = 0;
			int index = 0;
			do
			{
				if(start + index >= 0 && start + index < in.length()) // Do not check if the beginning of "search" is before the beginning of "in"
				{
					if(search.charAt(index)!= '*' && in.charAt(start + index) != search.charAt(index)) mismatch++;
				}
				index++;
			} while(index < search.length() && mismatch <= allowedMismatch);
			if(mismatch <= allowedMismatch)
			{
				if(start < res.offset) res.offset = start; // Search for bigger overlap
				res.overlap = (lE - lS);
				res.mismatch = mismatch;
			}
		}
		return res;
	}
	
	private static Result searchForPolyA(String in)
	{
		Result res = new Result();
		int maxPolyALength = 0;
		for(int i = 0; i < in.length(); i++) 
		{
			char c = in.charAt(i);
			if(c == 'A') res.polyALength++;
			else
			{
				if(maxPolyALength < res.polyALength) maxPolyALength = res.polyALength;
				res.polyALength = 0;
			}
		}
		res.polyALength = maxPolyALength;
		return res;
	}
	
	private static Result trimPolyA(String in)
	{
		Result res = new Result();
		res.offset = 0;
		int maxPolyALength = 0;
		for(int i = 0; i < in.length(); i++) 
		{
			char c = in.charAt(i);
			if(c == 'A') 
			{
				res.polyALength++;
				if(!res.isTrimmed && res.polyALength >= Parameters.polyALength) 
				{
					res.isTrimmed = true;
					if(maxPolyALength > res.polyALength) res.polyALength = maxPolyALength;
					res.outputRead = in.substring(0, res.offset);
					return res;
				}
			}
			else
			{
				if(maxPolyALength < res.polyALength) maxPolyALength = res.polyALength;
				res.offset = i;
				res.polyALength = 0;
			}
		}
		res.outputRead = in;
		res.polyALength = maxPolyALength;
		return res;
	}
	
	private static String autoFindBarcode(BufferedReader br) throws Exception
	{
		HashMap<String, Integer> barcodeOccurrence = new HashMap<String, Integer>();
		String line = br.readLine();
		long l = 1;
		while(line != null)
		{
			if(!line.startsWith("@")) 
			{
				br.close();
				System.err.println("Error while parsing FastQ file.");
				System.err.println("Line should start with a @, at l." + l);
				System.exit(-1);
			}
			String line2 = br.readLine().trim(); l++;
			br.readLine(); l++;
			br.readLine(); l++;
			String construct = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA***************" + "******" + "GT" + "ACTCTGCGTTGATACCACTGCTT" + "CCGCGGACAGGC" + "GTGTAGATCTCGGTGGTCGCCGTATCATT"; //everything after this should be discarded
			Result searchResult = searchFor(construct, line2, 10, 10);
			if(searchResult.mismatch != -1) // Barcode is from 45 to 50 
			{
				int barcodeStart = searchResult.offset + 45;
				if(barcodeStart >= 0 && barcodeStart < (line2.length() - 6))
				{
					String barcode = line2.substring(barcodeStart, barcodeStart + 6);
					if(!barcode.equals("AAAAAA") && !barcode.equals("GGGGGG"))
					{
						Integer count = barcodeOccurrence.get(barcode);
						if(count == null) count = 0;
						barcodeOccurrence.put(barcode, count + 1);
						if((count + 1) == 100) 
						{
							System.out.println("Barcode detected: " + reverseComplement(barcode) + " (>=100 times)");
							return barcode;
						}
					}
				}
			}
			line = br.readLine(); l++;
		}
		if(barcodeOccurrence.isEmpty()) return null;
		int max = 0;
		String maxbc = "";
		for(String barcode:barcodeOccurrence.keySet())
		{
			int count = barcodeOccurrence.get(barcode);
			if(count < max)
			{
				max = count;
				maxbc = barcode;
			}
		}
		System.out.println("The limit of 100 occurences was not reached.");
		System.out.println("Most occuring barcode is " + reverseComplement(maxbc) + " ("+max+" times)");
		return null;
	}
	
	private static BufferedReader readFileFromBeginning(String fastQFile) throws Exception
	{
		File fastQ = new File(fastQFile);
		String name = fastQ.getName();
		if(fastQ.getAbsolutePath().endsWith(".fastq"))
		{
			currentFilename = name.substring(0, name.lastIndexOf(".fastq"));
			return new BufferedReader(new FileReader(fastQ));
		}
		else if(fastQ.getAbsolutePath().endsWith(".fq")) 
		{
			currentFilename = name.substring(0, name.lastIndexOf(".fq"));
			return new BufferedReader(new FileReader(fastQ));
		}
		else if(fastQ.getAbsolutePath().endsWith(".fastq.gz"))
		{
			currentFilename = name.substring(0, name.lastIndexOf(".fastq.gz"));
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastQ))));
		}
		else if(fastQ.getAbsolutePath().endsWith(".fq.gz")) 
		{
			currentFilename = name.substring(0, name.lastIndexOf(".fq.gz"));
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastQ))));
		}
		else
		{
			System.err.println("The extension of the FastQ file is not recognized : " + fastQ.getAbsolutePath());
			System.err.println("It should be '.fastq', '.fq', '.fq.gz' or '.fastq.gz'");
			System.exit(-1);
		}
		return null;
	}

	private static String reverseComplement(String dna) 
	{
		String out = "";
		for(int i = dna.length() - 1; i >= 0; --i) 
		{
			char curr = dna.charAt(i);
			switch(curr)
			{
				case 'A': out += 'T'; break;
				case 'T': out += 'A'; break;
				case 'C': out += 'G'; break;
				case 'G': out += 'C'; break;
				case 'N': out += 'N'; break;
				default:
	                System.out.println("ERROR: Input is not a valid DNA Sequence.");
	                System.exit(-1);
			}
		}
		return out;
	}
	
	public static Excel trimFastQ(String fastQFile) throws Exception
	{
		long start = System.currentTimeMillis();
		Excel sheet = new Excel();
		
		BufferedReader br = readFileFromBeginning(fastQFile);
		GZIPOutputStream gzos = new GZIPOutputStream(new FileOutputStream(Parameters.outputFolder + currentFilename + ".trimmed.fastq.gz"));
		BufferedWriter bw_fq_trimmed = new BufferedWriter(new OutputStreamWriter(gzos, "UTF-8"));
		
		// Finding barcode for better construct trimming
		String barcode_RC = null;
		if(Parameters.autoFindBarcode)
		{
			barcode_RC = autoFindBarcode(br);
			br.close();
			if(barcode_RC == null) System.err.println("No barcode was detected during the automatic search... Switching to standard trimming mode.");
			br = readFileFromBeginning(fastQFile);
		}
		
		// Reading and trimming the fastq file
		String line = br.readLine();
		long l = 1;
		while(line != null)
		{
			if(!line.startsWith("@")) 
			{
				br.close();
				System.err.println("Error while parsing FastQ file " + fastQFile);
				System.err.println("Line should start with a @, at l." + l);
				System.exit(-1);
			}
			sheet.nbReads++;
			String line2 = br.readLine().trim(); l++;
			String line3 = br.readLine(); l++;
			String line4 = br.readLine(); l++;
			if(sheet.nbReads %Parameters.chunkSize == 0) System.out.println(sheet.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			
			Result searchPolyA = searchForPolyA(line2); // search for biggest polyA sequence (no allowed mismatch)
			polyADist_before[searchPolyA.polyALength]++;
			String construct = "******" + "GT" + "ACTCTGCGTTGATACCACTGCTT" + "CCGCGGACAGGC" + "GTGTAGATCTCGGTGGTCGCCGTATCATT"; //everything after this should be discarded
			if(barcode_RC != null) construct = barcode_RC + "GT" + "ACTCTGCGTTGATACCACTGCTT" + "CCGCGGACAGGC" + "GTGTAGATCTCGGTGGTCGCCGTATCATT"; //everything after this should be discarded

			Result searchResult = searchFor(construct, line2, 10, 10);
			if(searchResult.mismatch != -1)
			{
				sheet.nbContaminated++;
				contaminated[searchResult.overlap]++;
				if(searchResult.offset > 0)
				{
					String trimmedLine = line2.substring(0, searchResult.offset);
					Result trimPolyA = trimPolyA(trimmedLine);
					if(trimPolyA.isTrimmed) sheet.nbPolyATrimmed++;
					int le = trimPolyA.outputRead.length();
					if(le < Parameters.minReadLength) sheet.nbRemoved++;
					else
					{
						sheet.nbRemaininingReads++;
						polyADist_after[trimPolyA.polyALength]++;
						bw_fq_trimmed.write(line + "\n" + trimPolyA.outputRead + "\n" + line3 + "\n" + line4.substring(0, le) + "\n");
					}
				}
				else sheet.nbRemoved++;
			}
			else // No adapter trimming
			{
				Result trimPolyA = trimPolyA(line2);
				if(trimPolyA.isTrimmed) sheet.nbPolyATrimmed++;
				int le = trimPolyA.outputRead.length();
				if(le < Parameters.minReadLength) sheet.nbRemoved++;
				else
				{
					sheet.nbRemaininingReads++;
					polyADist_after[trimPolyA.polyALength]++;
					bw_fq_trimmed.write(line + "\n" + trimPolyA.outputRead + "\n" + line3 + "\n" + line4.substring(0, le) + "\n");
				}
			}
			line = br.readLine(); l++;
		}
		
		System.out.println(sheet.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		
		BufferedWriter bw_polyA = new BufferedWriter(new FileWriter(Parameters.outputFolder + currentFilename+ ".polyA.before.txt"));
		for(int i = 0; i < polyADist_before.length; i++) 
		{
			bw_polyA.write(i+"\t"+polyADist_before[i]+"\n");
			if(i >= Parameters.polyALength) sheet.nbPolyABefore += polyADist_before[i];
		}
		bw_polyA.close();
		
		bw_polyA = new BufferedWriter(new FileWriter(Parameters.outputFolder + currentFilename+ ".polyA.after.txt"));
		for(int i = 0; i < polyADist_after.length; i++) bw_polyA.write(i+"\t"+polyADist_after[i]+"\n");
		bw_polyA.close();
		
		BufferedWriter bw_contaminated = new BufferedWriter(new FileWriter(Parameters.outputFolder + currentFilename + ".adapterContamination.txt"));
		for(int i = 0; i < contaminated.length; i++) bw_contaminated.write(i+"\t"+contaminated[i]+"\n");
		bw_contaminated.close();
		
		bw_fq_trimmed.close();
		gzos.finish();
		gzos.close();
		
		System.out.println(sheet);
		return sheet;
	}
}

class Result
{
	int offset = Integer.MAX_VALUE;
	int overlap = 0;
	int mismatch = -1;
	String overlapString = "";
	boolean isTrimmed = false;
	int polyALength = 0;
	String outputRead = null;
}

