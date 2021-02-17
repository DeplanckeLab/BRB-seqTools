import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import model.Excel;
import model.GTF;
import model.Parameters;
import tools.Utils;

/**
 * @author Vincent Gardeux
 * @see vincent.gardeux@epfl.ch
 *
 */
public class BRBseqTools
{	
	// TODO Gérer les TMP files si meme folder
	// TODO Bug avec Tophat? Multiple Mapping?
	public static void main(String[] args) throws Exception
	{
		if(args.length < 1) Parameters.printHelp();
		else
		{
			String[] argsParsed = new String[args.length - 1];
			for(int i = 0; i < args.length - 1; i++) argsParsed[i] = args[i + 1];
			
			switch(args[0])
			{
				case "CreateDGEMatrix":
					System.out.println("BRBSeqTools " + Parameters.currentVersion + " [CreateDGEMatrix]\n");
					Parameters.loadDGE(argsParsed);
					Parameters.barcodes = Utils.readConfig();
					Parameters.BC1 = new ArrayList<String>();
					Parameters.barcodeIndex = new HashMap<String, Integer>();
					for(String B1:Parameters.barcodes) Parameters.BC1.add(B1);
					for(int i = 0; i < Parameters.BC1.size(); i++) Parameters.barcodeIndex.put(Parameters.BC1.get(i), i);
					Parameters.barcodeIndex.put("Unknown", Parameters.barcodeIndex.size());
					Utils.patterning();
					GTF.readGTF();
					Utils.readR1Fastq(false);
					DGEMatrixManager.readR2BAM();
					DGEMatrixManager.createOutputDGE();
					break;
				case "ExtractReadCountMatrix":
					System.out.println("BRBSeqTools " + Parameters.currentVersion + " [ExtractReadCountMatrix]\n");
					Parameters.loadExtractReadCountMatrix(argsParsed);
					Parameters.barcodes = Utils.readConfig();
					Parameters.BC1 = new ArrayList<String>();
					Parameters.barcodeIndex = new HashMap<String, Integer>();
					for(String B1:Parameters.barcodes) Parameters.BC1.add(B1);
					for(int i = 0; i < Parameters.BC1.size(); i++) Parameters.barcodeIndex.put(Parameters.BC1.get(i), i);
					System.out.println(Parameters.barcodeIndex.size() + " samples in barcode file.");
					Parameters.barcodeIndex.put("Unknown", Parameters.barcodeIndex.size());
					GTF.readGTF();
					ExtractReadCountMatrixManager.extract();
					break;
				case "AnnotateBAM":
					System.out.println("BRBSeqTools " + Parameters.currentVersion + " [AnnotateBAM]\n");
					Parameters.loadAnnoBAM(argsParsed);
					if(Parameters.inputConfigFile != null)
					{
						Parameters.barcodes = Utils.readConfig();
						Parameters.BC1 = new ArrayList<String>();
						Parameters.barcodeIndex = new HashMap<String, Integer>();
						for(String B1:Parameters.barcodes) Parameters.BC1.add(B1);
						for(int i = 0; i < Parameters.BC1.size(); i++) Parameters.barcodeIndex.put(Parameters.BC1.get(i), i);
						Parameters.barcodeIndex.put("Unknown", Parameters.barcodeIndex.size());
					}
					Utils.patterning();
					if(Parameters.inputGTFFile != null) GTF.readGTF();
					Utils.readR1Fastq(true);
					AnnotateBAMManager.readR2BAM();
					AnnotateBAMManager.annotate();
					break;
				case "Demultiplex":
					System.out.println("BRBSeqTools " + Parameters.currentVersion + " [Demultiplex]\n");
					Parameters.loadDemultiplex(argsParsed);		
					Parameters.barcodes = Utils.readConfig();
					Parameters.BC1 = new ArrayList<String>();
					Parameters.barcodeIndex = new HashMap<String, Integer>();
					for(String B1:Parameters.barcodes) Parameters.BC1.add(B1);
					for(int i = 0; i < Parameters.BC1.size(); i++) Parameters.barcodeIndex.put(Parameters.BC1.get(i), i);
					Parameters.barcodeIndex.put("Unknown", Parameters.barcodeIndex.size());
					Utils.patterning();
					DemultiplexingManager.demultiplex();
					break;
				case "Trim":
					System.out.println("BRBSeqTools " + Parameters.currentVersion + " [Trim]\n");
					Parameters.loadTrim(argsParsed);
					BufferedWriter bw_excel = new BufferedWriter(new FileWriter(Parameters.outputFolder + "brbseq.trimming.output.txt"));
					bw_excel.write("sample\tnbReads\tnbRemaininingReads\tnbContaminated\tnbPolyATrimmed\tnbRemoved\tnbPolyABefore\n");
					for(File fi:Parameters.fastqToTrim) 
					{
						System.out.println("\nTrimming " + fi.getName() + "...");
						long start = System.currentTimeMillis();
						Excel sheet = TrimmingManager.trimFastQ(fi.getAbsolutePath().replaceAll("\\\\", "/"));
						bw_excel.write(fi.getName());
						sheet.write(bw_excel);
						System.out.println("Trimmming done in " + Utils.toReadableTime(System.currentTimeMillis() - start));
					}
					bw_excel.close();
					break;
				default:
					System.err.println("This method is not implemented. Please use one of the following: [CreateDGEMatrix, AnnotateBAM, Demultiplex, Trim].");
					Parameters.printHelp();
			}
		}
	}
	

}

