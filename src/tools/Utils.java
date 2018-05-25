package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import model.GTF;
import model.Parameters;
import model.Read;
import model.TmpFileManager;

public class Utils 
{
	/**
	 * Read recursively a folder to find FastQ files
	 * @param directory root folder to search into
	 */
	public static void listFASTQfiles(File directory, ArrayList<File> files) 
	{
		File[] fList = directory.listFiles();
		for (File file : fList) 
		{
			if (file.isFile() && (file.getName().endsWith(".fq.gz") || file.getName().endsWith(".fastq.gz") || file.getName().endsWith(".fq") || file.getName().endsWith(".fastq"))) files.add(file);
			else if (file.isDirectory()) listFASTQfiles(file, files);
		}
	}
	
	/**
	 * Read FastQ file containing the mapping between read name, barcode and UMI.
	 * @param fastQ input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static BufferedReader readFastq(File fastQ) throws Exception
	{
		if(fastQ.getAbsolutePath().endsWith(".fastq") || fastQ.getAbsolutePath().endsWith(".fq"))
		{
			return new BufferedReader(new FileReader(fastQ));
		}
		else if(fastQ.getAbsolutePath().endsWith(".fastq.gz") || fastQ.getAbsolutePath().endsWith(".fq.gz"))
		{
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
	
	public static BufferedReader readGTF(File gtf) throws Exception
	{
		if(gtf.getAbsolutePath().endsWith(".gtf"))
		{
			return new BufferedReader(new FileReader(gtf));
		}
		else if(gtf.getAbsolutePath().endsWith(".gtf.gz"))
		{
			return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gtf))));
		}
		else
		{
			System.err.println("The extension of the GTF file is not recognized : " + gtf.getAbsolutePath());
			System.err.println("It should be '.gtf', or '.gtf.gz'");
			System.exit(-1);
		}
		return null;
	}
	
	/**
	 * Read FastQ file containing the mapping between read name, barcode, UMI and phreds
	 * @param fastQ input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static Read nextRead(BufferedReader br, boolean doQual) throws Exception
	{
		Read read = new Read();
		String header = br.readLine(); // Get first line = READNAME + INDEX
		if(header == null) return null;
		if(!header.startsWith("@")) throw new Exception("R1 fastq file has formatting issues");
		read.name = header.substring(1, header.indexOf(" "));
		String indexes = br.readLine().trim(); // Get second line = READ
		if(indexes.length() != Parameters.lengthBarcode) // Checking the length is OK
		{
			br.close();
			System.err.println("Error while parsing R1 FastQ file");
			System.err.println("Read found in FastQ has length " + indexes.length() + " while barcode pattern has length " + Parameters.lengthBarcode);
			System.exit(-1);
		}
		String barcode = null;
		try
		{
			if(Parameters.l1 != -1) barcode = indexes.substring(Parameters.barcode1Range[0], Parameters.barcode1Range[1]); // If there is a barcode to look for
			if(Parameters.UMILength != -1) read.UMI = indexes.substring(Parameters.UMIRange[0], Parameters.UMIRange[1]);
		}
		catch(IndexOutOfBoundsException boe)
		{
			br.close();
			System.err.println("Error while parsing R1 FastQ file");
			System.err.println("Index found in FastQ "+ indexes +" is not according to the pattern: " + Parameters.barcodePattern);
			System.exit(-1);
		}
		if(Parameters.l1 != -1) // If there is a barcode to look for
		{
			if(Parameters.inputConfigFile != null) // If there is a configuration file specified
			{
				// Check if barcode exists
				if(Parameters.BC1.contains(barcode)) { read.barcodeMatch = true; read.barcode = barcode;}
				else
				{
					String barcodeSaved = saveBarcode(barcode, Parameters.BC1, Parameters.nbAllowedDiff); // Saving
					if(barcodeSaved != null)
					{
						read.barcodeMatch = true;
						read.barcode = barcodeSaved;
					}
					else read.barcode = barcode; // Barcode not FOUND and not SAVED. By default barcodeSaved = false
				}
			}
			else read.barcode = barcode;
			if(Parameters.addReadGroup) // Do it in advance to directly create the correct header in output annotated BAM file (in case it is needed)
			{
				String[] tokens = read.name.split(":"); // TODO what if this does not work? Here the example that works is readName="NB500883:254:HL2K3BGX5:3:13509:21553:18111"
				String sampleName = "UNKNOWN";
				if(read.barcodeMatch) sampleName = Parameters.mappingBarcodeName.get(read.barcode).replaceAll("\\.", "_");
				String id = tokens[2] + "." + sampleName + "." + tokens[3];
				Parameters.readGroups.add(id);
			}
		}
		br.readLine(); // Third line = we don't care
		String line = br.readLine();  // Fourth line = qualities
		if(doQual)
		{
			line = line.trim();
			try
			{
				if(Parameters.l1 != -1) read.qualityBC = line.substring(Parameters.barcode1Range[0], Parameters.barcode1Range[1]);
				if(Parameters.UMILength != -1) read.qualityUMI = line.substring(Parameters.UMIRange[0], Parameters.UMIRange[1]);
			}
			catch(IndexOutOfBoundsException boe)
			{
				br.close();
				System.err.println("Error while parsing R1 FastQ file PHRED scores");
				System.exit(-1);
			}
		}
		return read;
	}
	
	/**
	 * Read R2 FastQ file containing the sequence read
	 * @param BufferedReader handle
	 * @return Read read
	 * @throws Exception
	 */
	public static Read nextDataRead(BufferedReader br) throws Exception
	{
		Read read = new Read();
		read.rawData = new String[4]; // 4 lines
		String header = br.readLine(); // Get first line = READNAME + INDEX
		if(header == null) return null;
		if(!header.startsWith("@")) throw new Exception("R2 fastq file has formatting issues");
		read.name = header.substring(1, header.indexOf(" "));
		read.rawData[0] = header; 
		read.rawData[1] = br.readLine(); // Get second line = READ
		read.rawData[2] = br.readLine(); // Get Third line = "+"
		read.rawData[3] = br.readLine(); // Get Fourth line = Quality
		return read;
	}
	
	public static String saveBarcode(String barcode, ArrayList<String> barcodes, int allowedDiff)
	{
		int[] diff = new int[barcodes.size()];
		for (int i = 0; i < barcode.length(); i++) 
		{
			for (int j = 0; j < barcodes.size(); j++) 
				if(barcode.charAt(i) != barcodes.get(j).charAt(i)) diff[j]++;
		}
		String sav = "";
		for (int i = 0; i < diff.length; i++) 
		{
			if(diff[i] == 0) 
			{
				sav = barcodes.get(i);
				break;
			}
			if(diff[i] <= allowedDiff)
			{
				if(!sav.equals(""))
				{
					//System.out.println(barcode + " cannot be saved because several barcodes can match it with " + allowedDiff + " errors.");
					sav = "";
					break;
				}
				sav = barcodes.get(i);
			}
		}
		if(sav.equals("")) return null;
		return sav;
	}
	
	public static HashSet<String> getOverlappingGenes(String chr, int start, int end, Cigar c, boolean readNegativeStrandFlag) throws Exception
	{
		HashSet<String> res = new HashSet<>();
		List<CigarElement> l = c.getCigarElements();
		int s = start;
		for(CigarElement cigar:l)
		{
			switch(cigar.getOperator())
			{
				case M:
					res.addAll(GTF.findOverlappingGenes(chr, s, s + cigar.getLength() - 1, readNegativeStrandFlag)); // -1 Because the last letter is at the index before
					s += cigar.getLength();
					break;
				case N:
					s += cigar.getLength();
					break;
				case D:
					s += cigar.getLength();
					break;
				case EQ:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case H:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case I:
					// Do nothing
					break;
				case P:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
				case S:
					// Do nothing (alignment start is after any S & alignment end before any S)
					break;
				case X:
					System.out.println("CIGAR = " + c);
					return new HashSet<>(); // TODO
			}
		}
		
		if(s != end + 1) throw new Exception("Error while reading CIGAR");
		return res;
	}
	
	/**
	 *  Reading reads barcodes/UMI from the R1 fastq file to map the UMI/barcode with the read name (lost after alignment)
	 * @throws Exception Yes I know...
	 */
	public static void readR1Fastq(boolean readQual) throws Exception
	{
		System.out.println("\nReading reads barcodes/UMI from the R1 fastq file...");
		BufferedReader br = Utils.readFastq(Parameters.inputFastQFileR1);
		TreeSet<String> lines = new TreeSet<String>();
		int noBarcodeMatch = 0;
		Long start = System.currentTimeMillis();
		Parameters.nbTmpFastq = 1; // Number of temporary files
		Read read = Utils.nextRead(br, readQual);
		while(read != null)
		{
			Parameters.nbReads++;
			if(Parameters.nbReads %Parameters.chunkSize == 0) 
			{
				TmpFileManager.createTmpFile(Parameters.tmpFolder + "fastq."+Parameters.nbTmpFastq+".tmp", lines);
				Parameters.nbTmpFastq++;
				lines = new TreeSet<String>();
				System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			if(!read.barcodeMatch && Parameters.inputConfigFile != null) { noBarcodeMatch++; read.barcode = "Unknown";}
			if(!readQual) lines.add(read.name + "\t" + read.barcode + "\t" + read.UMI);
			else lines.add(read.name + "\t" + read.barcode + "\t" + read.UMI + "\t" + read.qualityBC + "\t" + read.qualityUMI);
			read = Utils.nextRead(br, readQual);
		}
		br.close();
				
		TmpFileManager.createTmpFile(Parameters.tmpFolder + "fastq."+Parameters.nbTmpFastq+".tmp", lines);
		System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");		
		System.out.println("Created " + Parameters.nbTmpFastq + " temporary fastq files");
		
		if(Parameters.inputConfigFile != null) System.out.println(noBarcodeMatch + " reads have no matching barcodes (" + Parameters.pcFormatter.format(((float)noBarcodeMatch / Parameters.nbReads) * 100) + "%)");
	}
	
	/**
	 * 
	 * @param configFile Required to map the barcode with the sample name
	 * @return Unique set of all barcodes
	 * @throws Exception Yes I know...
	 */
	public static HashSet<String> readConfig() throws Exception
	{
		HashSet<String> res = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(Parameters.inputConfigFile));
		String line = br.readLine();
		int indexes[] = new int[3];
		Arrays.fill(indexes, -1);
		String[] header = line.split("\t");
		for (int i = 0; i < header.length; i++) 
		{
			if(header[i].trim().equals("B1")) indexes[0] = i;
			if(header[i].trim().equals("Name")) { indexes[1] = i; Parameters.mappingBarcodeName = new HashMap<String, String> (); }
		}
		if(indexes[0] == -1)
		{
			System.err.println("B1 column not found in Config file (Required).");
			System.exit(-1);
		}
		if(indexes[1] == -1)
		{
			System.err.println("Name column not found in Config file (Required).");
			System.exit(-1);
		}
		line = br.readLine();
		while(line != null)
		{
			String[] barcodes = line.split("\t");
			if(barcodes.length >= 1)
			{
				String B1 = barcodes[indexes[0]].trim();
				if(indexes[1] != -1) Parameters.mappingBarcodeName.put(B1, barcodes[indexes[1]].trim());
				res.add(B1);
			}
			line = br.readLine();
		}
		br.close();
		// Quality Check
		Parameters.l1 = -1;
		System.out.print("Config: B1 contains " + res.size() + " barcodes: [ ");
		for(String B1:res)
		{
			if(Parameters.l1 == -1) Parameters.l1 = B1.length();
			if(Parameters.l1 != B1.length())
			{
				System.err.println("Error: All the 'B1' barcodes should have the same length. Check your config file.");
				System.exit(-1);
			}
			System.out.print(B1 + " ");
		}
		System.out.println("]\nConfig: B1 length [#nucleotides] = " + Parameters.l1);
		return res;
	}
	
	/**
	 * Parsing the barcode pattern and check validity
	 * @throws Exception Yes I know...
	 */
	public static void patterning() throws Exception
	{
		Parameters.lengthBarcode = 0;
		System.out.println("\nAnalyzing barcode Pattern...");
		for(int j = 0; j < Parameters.barcodePattern.length(); j++) 
		{
			char c = Parameters.barcodePattern.charAt(j);
			switch(c)
			{
				case 'B':
					Parameters.barcode1Range[0] = Parameters.lengthBarcode;
					Parameters.lengthBarcode += Parameters.l1;
					Parameters.barcode1Range[1] = Parameters.lengthBarcode;
					break;					
				case 'U':
					Parameters.UMIRange[0] = Parameters.lengthBarcode;
					Parameters.lengthBarcode += Parameters.UMILength;
					Parameters.UMIRange[1] = Parameters.lengthBarcode;
					break;
				case '?':
					Parameters.lengthBarcode++;
					break;
				default:
					System.err.println(c+" does not correspond to any authorized pattern character (?, B, U). Aborted.");
					System.exit(-1);
			}
		}
		System.out.println("According to barcode pattern, reads of R1 FastQ file should contain "+Parameters.lengthBarcode+" characters.");
	}
	
	public static String[] sortKeys(Map<String, Integer> map)
	{
		String[] keys = map.keySet().toArray(new String[map.keySet().size()]);
		Arrays.sort(keys);
		return keys;
	}
	
	public static String toReadableTime(long ms)
	{
		if(ms < 1000) return ""+ms+" ms";
		long s = ms / 1000;
		ms = ms % 1000;
		if(s < 60) return s+" s "+ms+" ms";
		long mn = s / 60;
		s = s % 60;
		if(mn < 60) return mn +" mn "+s+" s "+ms+" ms";
		long h = mn / 60;
		mn = mn % 60;
		if(h < 24) return h +" h "+ mn +" mn "+s+" s "+ms+" ms";
		long d = h / 24;
		h = h % 24;
		return d+ " d " + h +" h "+ mn +" mn "+s+" s "+ms+" ms";
	}
}
