import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import model.Parameters;
import model.Read;
import tools.Utils;

public class DemultiplexingManager 
{
	public static HashMap<String, BufferedWriter> bw = null; // handle on all fastq files to generate
	
	public static void demultiplex() throws Exception
	{
		Long start = System.currentTimeMillis();
		System.out.println("\nStarting demultiplexing...");
		Parameters.nbReads = 0;
		int noBarcodeMatch = 0;
		HashMap<String, Integer> nbReadsPerBarcode = new HashMap<>();
		for (String B1:Parameters.BC1) nbReadsPerBarcode.put(B1, 0);
		
		// Create Writing Buffers
		bw = new HashMap<String, BufferedWriter>();
		for (String barcode1:Parameters.barcodes) bw.put(barcode1, new BufferedWriter(new OutputStreamWriter( new GZIPOutputStream(new FileOutputStream(Parameters.outputFolder + Parameters.mappingBarcodeName.get(barcode1) + ".fastq.gz")), "UTF-8")));
		bw.put("undetermined", new BufferedWriter(new OutputStreamWriter( new GZIPOutputStream(new FileOutputStream(Parameters.outputFolder + "undetermined.fastq.gz")), "UTF-8")));
		//UndefinedFiles.add(Parameters.outputFolder + "undefinedBarcodes.txt");
		bw.put("stats", new BufferedWriter(new FileWriter(Parameters.outputFolder + "stats.txt")));

		// Read R1 & R2 files
		BufferedReader br = Utils.readFastq(Parameters.inputFastQFileR1);
		BufferedReader br2 = Utils.readFastq(Parameters.inputFastQFileR2);
		
		// Initialize by reading the 2 first reads of both R1 & R2 files
		Read read1 = Utils.nextRead(br, false);
		Read read2 = Utils.nextDataRead(br2);
		
		// Check the separator for writing UMI in output fastq files
		char lastChar = read2.rawData[0].charAt(read2.rawData[0].length() - 1);
		if(lastChar == 'G' || lastChar == 'A' || lastChar == 'C' || lastChar == 'T') Parameters.separator = "+";
		else if(lastChar == ':') Parameters.separator = "";
		else if(lastChar == 'N')			
		{
			lastChar = read2.rawData[0].charAt(read2.rawData[0].length() - 2); // check character before
			if(lastChar != ':') Parameters.separator = "+";
		}
		
		// Read the 2 fastq files
		while(read1 != null) // For each read
		{
			Parameters.nbReads++;
			if(Parameters.nbReads %Parameters.chunkSize == 0) System.out.println(Parameters.nbReads + " reads were processed from R1 & R2 fastq files [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			
			if(!read1.barcodeMatch) 
			{
				noBarcodeMatch++; 
				read2.write(bw.get("undetermined"), read1.UMI);
			}
			else
			{
				nbReadsPerBarcode.put(read1.barcode, nbReadsPerBarcode.get(read1.barcode) + 1);
				if(read1.name.equals(read2.name)) read2.write(bw.get(read1.barcode), read1.UMI);
				else
				{
					System.err.println("ERROR: Read names do not match: R1= " + read1.name + " & R2= " + read2.name + "\nOf note: Both fastq (R1 & R2) should contain reads sorted in the same order (should be the case by default)");
					System.exit(-1);
				}
			}
			read1 = Utils.nextRead(br, false);
			read2 = Utils.nextDataRead(br2);
		}
	    br.close();
	    br2.close();
		System.out.println(Parameters.nbReads + " reads were processed from R1 & R2 fastq files [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");

		// Summary Files
		bw.get("stats").write("Demultiplexing Parameters\n");
		bw.get("stats").write("R1 fastq file\t" + Parameters.inputFastQFileR1 + "\n");
		bw.get("stats").write("R2 fastq file\t" + Parameters.inputFastQFileR2 + "\n");
		bw.get("stats").write("Path of barcode/samplename mapping file\t" + Parameters.inputConfigFile + "\n");
		bw.get("stats").write("Nb allowed barcode mismatch\t" + Parameters.nbAllowedDiff + "\n");
		bw.get("stats").write("Output folder\t" + Parameters.outputFolder + "\n");
		bw.get("stats").write("Barcode pattern in R1\t" + Parameters.barcodePattern + "\n");
		if(Parameters.UMILength == -1) bw.get("stats").write("UMI length\tNO_UMI\n");
		else bw.get("stats").write("UMI length\t" + Parameters.UMILength + "\n");
		bw.get("stats").write("\nDemultiplexing Summary\n");
		bw.get("stats").write("Summary\tNb_Reads\tPercent\n");
		bw.get("stats").write("Total number of reads\t" + Parameters.nbReads + "\t100\n");
		bw.get("stats").write("Undetermined Barcodes\t" + noBarcodeMatch + "\t" + Parameters.pcFormatter.format(((double)noBarcodeMatch / Parameters.nbReads) * 100) + "\n");
		bw.get("stats").write("Determined Barcodes\t" + (Parameters.nbReads - noBarcodeMatch) + "\t" + Parameters.pcFormatter.format(((double)(Parameters.nbReads - noBarcodeMatch) / Parameters.nbReads) * 100) + "\n");
		bw.get("stats").write("\nDemultiplexing Details\n");
		bw.get("stats").write("Barcode\tDemult_Reads\tPercent\n");
		for (String B1:Parameters.BC1) bw.get("stats").write(Parameters.mappingBarcodeName.get(B1) + "\t" + nbReadsPerBarcode.get(B1) + "\t" + Parameters.pcFormatter.format(((double)nbReadsPerBarcode.get(B1) / Parameters.nbReads) * 100) + "\n");

		// Close Writing Buffers
		for (String key:bw.keySet()) bw.get(key).close();
	}
}
