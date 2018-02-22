package model;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import tools.Utils;

enum Strand{NO, YES, REVERSE};

public class Parameters 
{
	// Input parameters
	public static String tmpFolder = null;
	public static String outputFolder = null;
	public static File inputConfigFile = null;
	public static File inputFastQFileR1 = null;
	public static File inputFastQFileR2 = null;
	public static File inputBAMFileR2 = null;
	public static File inputGTFFile = null;
	public static ArrayList<File> fastqToTrim = null;
	public static String barcodePattern = "BU";
	public static int nbAllowedDiff = 1;
	public static int UMILength = -1;
	public static int polyALength = 6;
	public static int minReadLength = 10;
	public static Strand stranded = Strand.YES; 
	public static long chunkSize = 1000000;
	public static boolean autoFindBarcode = false;
	public static String separator = ":";
	
	// Computed variables
	public static int noFeature = 0;
	public static int notUnique = 0;
	public static int ambiguous = 0;
	public static int mapped = 0;
	public static int unmapped = 0;
	public static int toolowAqual = 0;
	public static double readLength = 0;
	public static int nbTmpFastq = 0;
	public static int nbTmpBAM = 0;
	public static long nbReads = 0;
	public static int[] barcode1Range = new int[2];
	public static int[] UMIRange = new int[2];
	public static DecimalFormat pcFormatter = new DecimalFormat("##.##");
	public static HashSet<String> barcodes = null;
	public static ArrayList<String> BC1;
	public static int lengthBarcode = 0;
	public static HashMap<String, String> mappingBarcodeName = null;
	public static HashMap<String, String> mappingGeneIdGeneName = null;
	public static int l1 = -1;
	public static HashMap<String, Integer> geneIndex = null;
	public static HashMap<String, Integer> barcodeIndex = null;
	
	public static void loadDemultiplex(String[] args) throws Exception
	{
		boolean barcodeSpecified = false;
		if(args.length == 0)
		{
			printHelpDemultiplex();
			System.exit(0);
		}
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) throw new Exception(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-c":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputConfigFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-c' option should be followed by one config file path. " + e.getMessage() + ".");
							System.exit(-1);
						}
						break;
					case "-n":
						i++;
						try
						{
							nbAllowedDiff = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-n' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-UMI":
						i++;
						try
						{
							UMILength = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-UMI' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-r1":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputFastQFileR1 = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-r1' option should be followed by R1 FastQ file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-r2":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputFastQFileR2 = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-r2' option should be followed by R2 FastQ file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-p":
						i++;
						try
						{
							barcodeSpecified = true;
							barcodePattern = args[i];
							if(!barcodePattern.contains("B")) throw new Exception("'B' not found in R1 fastq Pattern (i.e. position of barcode)");
						}
						catch(Exception e)
						{
							System.err.println("The '-p' option should be followed by a valid barcode pattern. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(!barcodeSpecified) System.out.println("No barcode pattern specified (with '-p' option). Barcode pattern is set as default: '"+barcodePattern+"'");
		if(inputFastQFileR1 == null)
		{
			System.err.println("Please use '-r1' option to specify R1 FastQ file with barcode/UMI information");
			System.exit(-1);
		}
		if(inputFastQFileR2 == null)
		{
			System.err.println("Please use '-r2' option to specify R2 FastQ file with read sequence information");
			System.exit(-1);
		}
		if(inputConfigFile == null)
		{
			System.err.println("Please use '-c' option to specify a config file.");
			System.exit(-1);
		}
		if(UMILength == -1 && barcodePattern.contains("U"))
		{
			System.err.println("Your barcode pattern contains UMI, you should specify its length with the -UMI option.");
			System.exit(-1);
		}
		if(UMILength != -1 && !barcodePattern.contains("U"))
		{
			System.err.println("You specified a UMI length but your barcode pattern does not contain 'U', you should specify where to find the UMI in R1");
			System.exit(-1);
		}
		System.out.println("Config: Barcode Pattern = " + barcodePattern);
		System.out.println("Config: Nb Allowed Diff = " + nbAllowedDiff);
		if(outputFolder == null)
		{
			String path = inputConfigFile.getAbsolutePath();
			path = path.replaceAll("\\\\", "/");
			path = path.substring(0, path.lastIndexOf("/"));
			System.out.println("No output Folder is specified, using default one: \""+path+"\". You can specify an output path by using the '-o' option.");
			outputFolder = path;
		}
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
	}
	
	public static void loadDGE(String[] args) throws Exception
	{
		boolean barcodeSpecified = false;
		if(args.length == 0)
		{
			printHelpDGE();
			System.exit(0);
		}
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) throw new Exception(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-c":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputConfigFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-c' option should be followed by one config file path. " + e.getMessage() + ".");
							System.exit(-1);
						}
						break;
					case "-s":
						i++;
						switch(args[i])
						{
						case "no":
							Parameters.stranded = Strand.NO;
							break;
						case "yes":
							Parameters.stranded = Strand.YES;
							break;
						case "reverse":
							Parameters.stranded = Strand.REVERSE;
							break;
						default:
							System.err.println("The '-s' option should be followed by one of the following parameters: [no, yes, reverse].");
							System.exit(-1);
						}
						break;
					case "-n":
						i++;
						try
						{
							nbAllowedDiff = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-n' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-UMI":
						i++;
						try
						{
							UMILength = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-UMI' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-f":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputFastQFileR1 = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-f' option should be followed by R1 FastQ file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-b":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputBAMFileR2 = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-b' option should be followed by aligned BAM file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-p":
						i++;
						try
						{
							barcodeSpecified = true;
							barcodePattern = args[i];
							if(!barcodePattern.contains("B")) throw new Exception("B not found in R1 fastq Pattern (i.e. position of barcode)");
						}
						catch(Exception e)
						{
							System.err.println("The '-p' option should be followed by a valid barcode pattern. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-chunkSize":
						i++;
						try
						{
							chunkSize = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-chunkSize' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-gtf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputGTFFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-gtf' option should be followed by GTF file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-t":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No Directory at path " + args[i]);
							if(!c.isDirectory()) throw new Exception(args[i] + " is not a directory");
							tmpFolder = c.getAbsolutePath();
						}
						catch(Exception e)
						{
							System.err.println("The '-t' option should be followed by the path of a folder to store temporary files. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(!barcodeSpecified) System.out.println("No barcode pattern specified (with '-p' option). Barcode pattern is set as default: '"+barcodePattern+"'");
		if(inputFastQFileR1 == null)
		{
			System.err.println("Please use '-f' option to specify R1 FastQ file with barcode/UMI information");
			System.exit(-1);
		}
		if(inputBAMFileR2 == null)
		{
			System.err.println("Please use '-b' option to specify the path of the aligned BAM file");
			System.exit(-1);
		}
		if(inputConfigFile == null)
		{
			System.err.println("Please use '-c' option to specify a config file.");
			System.exit(-1);
		}
		if(inputGTFFile == null)
		{
			System.err.println("Please use '-gtf' option to specify the path of the GTF file");
			System.exit(-1);
		}
		if(UMILength == -1 && barcodePattern.contains("U"))
		{
			System.err.println("Your barcode pattern contains UMI, you should specify its length with the -UMI option.");
			System.exit(-1);
		}
		if(UMILength != -1 && !barcodePattern.contains("U"))
		{
			System.err.println("You specified a UMI length but your barcode pattern does not contain 'U', you should specify where to find the UMI in R1");
			System.exit(-1);
		}
		System.out.println("Barcode Pattern = " + barcodePattern);
		System.out.println("Nb Allowed Diff = " + nbAllowedDiff);
		System.out.println("ChunkSize = " + chunkSize + " i.e. no more than " + chunkSize + " reads will be stored in RAM.");
		System.out.println("Stranded = " + Parameters.stranded);
		if(outputFolder == null)
		{
			String path = inputConfigFile.getAbsolutePath();
			path = path.replaceAll("\\\\", "/");
			path = path.substring(0, path.lastIndexOf("/"));
			System.out.println("No output Folder is specified, using default one: \""+path+"\". You can specify an output path by using the '-o' option.");
			outputFolder = path;
		}
		if(tmpFolder == null) 
		{
			System.out.println("No folder for temporary generated files is specified, using output folder as default: \""+outputFolder+"\". You can specify a temporary folder by using the '-t' option.");
			tmpFolder = outputFolder;
		}
		tmpFolder = tmpFolder.replaceAll("\\\\", "/");
		if(!tmpFolder.endsWith("/")) tmpFolder+= "/";
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
	}
	
	public static void loadAnnoBAM(String[] args) throws Exception
	{
		boolean strandSpecified = false;
		boolean mismatchSpecified = false;
		boolean barcodeSpecified = false;
		if(args.length == 0)
		{
			printHelpAnnoBAM();
			System.exit(0);
		}
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) throw new Exception(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-UMI":
						i++;
						try
						{
							UMILength = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-UMI' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-c":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputConfigFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-c' option should be followed by one config file path. " + e.getMessage() + ".");
							System.exit(-1);
						}
						break;
					case "-BC":
						i++;
						try
						{
							l1 = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-BC' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-s":
						i++;
						strandSpecified = true;
						switch(args[i])
						{
						case "no":
							Parameters.stranded = Strand.NO;
							break;
						case "yes":
							Parameters.stranded = Strand.YES;
							break;
						case "reverse":
							Parameters.stranded = Strand.REVERSE;
							break;
						default:
							System.err.println("The '-s' option should be followed by one of the following parameters: [no, yes, reverse].");
							System.exit(-1);
						}
						break;
					case "-n":
						i++;
						mismatchSpecified = true;
						try
						{
							nbAllowedDiff = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-n' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-f":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputFastQFileR1 = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-f' option should be followed by R1 FastQ file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-b":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputBAMFileR2 = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-b' option should be followed by aligned BAM file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-p":
						i++;
						try
						{
							barcodeSpecified = true;
							barcodePattern = args[i];
						}
						catch(Exception e)
						{
							System.err.println("The '-p' option should be followed by a valid barcode pattern. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-chunkSize":
						i++;
						try
						{
							chunkSize = Long.parseLong(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-chunkSize' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-gtf":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No file at path " + args[i]);
							if(!c.isFile()) throw new Exception(args[i] + " is not a file");
							inputGTFFile = c;
						}
						catch(Exception e)
						{
							System.err.println("The '-gtf' option should be followed by GTF file path. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-t":
						i++;
						try
						{
							File c = new File(args[i]);
							if(!c.exists()) throw new Exception("No Directory at path " + args[i]);
							if(!c.isDirectory()) throw new Exception(args[i] + " is not a directory");
							tmpFolder = c.getAbsolutePath();
						}
						catch(Exception e)
						{
							System.err.println("The '-t' option should be followed by the path of a folder to store temporary files. " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(!barcodeSpecified) System.out.println("No barcode pattern specified (with '-p' option). Barcode pattern is set as default: '"+barcodePattern+"'");
		if(inputFastQFileR1 == null)
		{
			System.err.println("Please use '-f' option to specify R1 FastQ file with barcode/UMI information");
			System.exit(-1);
		}
		if(inputBAMFileR2 == null)
		{
			System.err.println("Please use '-b' option to specify the path of the aligned BAM file");
			System.exit(-1);
		}
		if(inputConfigFile == null) 
		{
			System.out.println("No barcode mapping file specified, barcodes will be annotated as in R1 (no mismatch correction)");
			if(mismatchSpecified)
			{
				System.err.println("You used the '-n' option for mapping barcodes (with mismatch correction). But you did not input a barcode mapping file. Whether remove '-n' option, or add a barcode mapping file with '-c'.");
				System.exit(-1);
			}
		}
		else
		{
			System.out.println("A barcode mapping file was specified, barcodes will be annotated as in R1 (with mismatch correction) or 'Unknown' if not mapping to any barcode.");
			System.out.println("Nb Allowed Diff = " + nbAllowedDiff);
		}
		if(barcodePattern.contains("B"))
		{
			if(l1 == -1 && inputConfigFile == null)
			{
				System.err.println("You did not specify a config file, nor a barcode length, but the barcode pattern contains 'B'. Please use one of these: '-BC' option to specify a barcode length, or '-c' option to specify a barcode mapping file");
				System.exit(-1);
			}
		}
		else 
		{
			if(l1 != -1) { System.err.println("Your barcode pattern does not contain any 'B' but you used the '-BC' option to specify a barcode length. Remove the '-BC' option, or change your barcorde pattern."); System.exit(-1);}
			if(inputConfigFile != null) { System.err.println("Your barcode pattern does not contain any 'B' but you used the '-c' option to specify a barcode mapping file. Remove the '-c' option, or change your barcorde pattern."); System.exit(-1);}
		}
		if(inputGTFFile == null)
		{
			System.out.println("No GTF file specified using '-gtf' option. Output annotated BAM will not contain annotated gene");
			if(strandSpecified) 
			{
				System.err.println("You used the '-s' option for annotating genes only on same strands. But you did not input a GTF file. Whether remove '-s' option, or add a GTF file with '-gtf'.");
				System.exit(-1);
			}
		}
		else
		{
			System.out.println("A GTF file was specified using the '-gtf' option. Output annotated BAM will contain annotated gene with 'CO' TAG");
			System.out.println("Stranded = " + Parameters.stranded);
		}
		if(UMILength == -1 && barcodePattern.contains("U"))
		{
			System.err.println("Your barcode pattern contains UMI, you should specify its length with the -UMI option.");
			System.exit(-1);
		}
		if(UMILength != -1 && !barcodePattern.contains("U"))
		{
			System.err.println("You specified a UMI length but your barcode pattern does not contain 'U', you should specify where to find the UMI in R1");
			System.exit(-1);
		}
		System.out.println("Barcode Pattern = " + barcodePattern);
		System.out.println("ChunkSize = " + chunkSize + " i.e. no more than " + chunkSize + " reads will be stored in RAM.");
		if(outputFolder == null)
		{
			String path = inputConfigFile.getAbsolutePath();
			path = path.replaceAll("\\\\", "/");
			path = path.substring(0, path.lastIndexOf("/"));
			System.out.println("No output Folder is specified, using default one: \""+path+"\". You can specify an output path by using the '-o' option.");
			outputFolder = path;
		}
		if(tmpFolder == null) 
		{
			System.out.println("No folder for temporary generated files is specified, using output folder as default: \""+outputFolder+"\". You can specify a temporary folder by using the '-t' option.");
			tmpFolder = outputFolder;
		}
		tmpFolder = tmpFolder.replaceAll("\\\\", "/");
		if(!tmpFolder.endsWith("/")) tmpFolder+= "/";
		outputFolder = outputFolder.replaceAll("\\\\", "/");
		if(!outputFolder.endsWith("/")) outputFolder+= "/";
		new File(outputFolder).mkdirs();
	}
	
	public static void loadTrim(String[] args) throws Exception
	{
		if(args.length == 0)
		{
			printHelpTrim();
			System.exit(0);
		}
		for(int i = 0; i < args.length; i++) 
		{
			if(args[i].startsWith("-"))
			{
				switch(args[i])
				{
					case "-o":
						i++;
						outputFolder = args[i];
						outputFolder = outputFolder.replaceAll("\\\\", "/");
						File f = new File(outputFolder);
						if(!f.exists()) 
						{
							System.out.println("Output folder does not exist. Creating it");
							f.mkdirs();
						}
						else if(!f.isDirectory()) throw new Exception(outputFolder + " is not a folder.");
						if(!outputFolder.endsWith("/")) outputFolder+= "/";
						break;
					case "-uniqueBarcode":
						autoFindBarcode = true;
						break;
					case "-polyA":
						i++;
						try
						{
							polyALength = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-polyA' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-minLength":
						i++;
						try
						{
							minReadLength = Integer.parseInt(args[i]);
						}
						catch(NumberFormatException nfe)
						{
							System.err.println("The '-minLength' option should be followed by an Integer. You entered " + args[i]);
							System.exit(-1);
						}
						break;
					case "-f":
						i++;
						try
						{
							String filename = args[i];
							filename = filename.replaceAll("\\\\", "/");
							f = new File(filename);
							if(!f.exists()) throw new Exception("No file at path " + filename);
							fastqToTrim = new ArrayList<File>();
							if(f.isDirectory()) 
							{
								Utils.listFASTQfiles(f, fastqToTrim);
								System.out.println(fastqToTrim.size() + " FASTQ files were found recursively at path: " + filename);
								tmpFolder = f.getAbsolutePath(); // in case no output folder is specified, I keep the containing folder as outputfolder
								tmpFolder = tmpFolder.replaceAll("\\\\", "/");
								if(!tmpFolder.endsWith("/")) tmpFolder+= "/";
							}
							else 
							{
								fastqToTrim.add(f);
								tmpFolder = f.getAbsolutePath(); // in case no output folder is specified, I keep the file location as outputfolder
								tmpFolder = tmpFolder.replaceAll("\\\\", "/");
								tmpFolder = tmpFolder.substring(0, tmpFolder.lastIndexOf("/"));
								if(!tmpFolder.endsWith("/")) tmpFolder+= "/";
							}

						}
						catch(Exception e)
						{
							System.err.println("The '-f' option should be followed by a FastQ file path (or containing folder path). " + e.getMessage() + ". You entered " + args[i]);
							System.exit(-1);
						}
						break;
					default:
						System.err.println("Unused argument: " + args[i]);
				}
			}
		}
		if(fastqToTrim == null || fastqToTrim.isEmpty())
		{
			System.err.println("Please use '-f' option to specify the path of a FastQ file (or containing folder for multiple fastq processing)");
			System.exit(-1);
		}
		if(outputFolder == null)
		{
			System.out.println("No output Folder is specified, using default one: \""+tmpFolder+"\".\nYou can specify an output path by using the '-o' option.");
			outputFolder = tmpFolder;
		}
	}
	
	private static void printHelpDGE()
	{
		System.out.println("BRBseq-tools 1.1 [CreateDGEMatrix]\n\nOptions:");
		System.out.println("-f %s \t\tPath of R1 FastQ file.");
		System.out.println("-b %s \t\tPath of R2 aligned BAM file [do not need to be sorted or indexed].");
		System.out.println("-c %s \t\tPath of Barcode/Samplename mapping file.");
		System.out.println("-gtf %s \tPath of GTF file.");
		System.out.println("-s %s \t\tDo you want to count only reads falling on same strand than gene? [no, yes, reverse] (default = yes since BRB-seq is stranded protocol).");
		System.out.println("-n %i \t\tNumber of allowed difference with the barcode [ambiguous reads will be automatically discarded].");
		System.out.println("-o %s \t\tOutput folder");
		System.out.println("-t %s \t\tPath of existing folder for storing temporary files");
		System.out.println("-chunkSize %i\tMaximum number of reads to be stored in RAM (default = 1000000)");
		System.out.println("-p %s \t\tBarcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).\n\t\t\t'B' [required] is used for specifying the barcode position.\n\t\t\t'U' [optional] can be used for specifying a UMI value position.\n\t\t\t'?' [optional] can be used to ignore specific nucleotides.");
		System.out.println("-UMI %i \tIf your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI.");
	}
	
	private static void printHelpAnnoBAM()
	{
		System.out.println("BRBseq-tools 1.1 [AnnotateBAM]\n\nOptions:");
		System.out.println("-f %s \t\tPath of R1 FastQ file.");
		System.out.println("-b %s \t\tPath of R2 aligned BAM file [do not need to be sorted or indexed].");
		System.out.println("-c %s \t\tPath of Barcode/Samplename mapping file.");
		System.out.println("-gtf %s \tPath of GTF file (for UMI counting).");
		System.out.println("-s %s \t\tDo you want to annotate genes only in case where reads are falling on same strand? [no, yes, reverse] (default = yes since BRB-seq is stranded protocol).");
		System.out.println("-n %i \t\tNumber of allowed difference with the barcode.");
		System.out.println("-o %s \t\tOutput folder");
		System.out.println("-t %s \t\tPath of existing folder for storing temporary files");
		System.out.println("-chunkSize %i\tMaximum number of reads to be stored in RAM (default = 1000000)");
		System.out.println("-p %s \t\tBarcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).\n\t\t\t'B' [required] is used for specifying the barcode position.\n\t\t\t'U' [optional] can be used for specifying a UMI value position.\n\t\t\t'?' [optional] can be used to ignore specific nucleotides.");
		System.out.println("-UMI %i \tYou should specify this parameter as the length of the UMI (required, since the duplicates are removed thanks to UMIs).");
	}
	
	private static void printHelpDemultiplex()
	{
		System.out.println("BRBseq-tools 1.1 [Demultiplex]\n\nOptions:");
		System.out.println("-r1 %s \t\tPath of R1 FastQ files (containing barcode and optionally UMIs).");
		System.out.println("-r2 %s \t\tPath of R2 FastQ files (containing read sequence).");
		System.out.println("-c %s \t\tPath of Barcode/Samplename mapping file.");
		System.out.println("-n %i \t\tNumber of allowed difference with the barcode [ambiguous reads will be automatically discarded].");
		System.out.println("-o %s \t\tOutput folder");
		System.out.println("-p %s \t\tBarcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).\n\t\t\t'B' [required] is used for specifying the barcode position.\n\t\t\t'U' [optional] can be used for specifying a UMI value position.\n\t\t\t'?' [optional] can be used to ignore specific nucleotides.");
		System.out.println("-UMI %i \tIf your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI.");
	}
	
	private static void printHelpTrim()
	{
		System.out.println("BRBseq-tools 1.1 [Trim]\n\nOptions:");
		System.out.println("-f %s \t\tPath of FastQ file to trim (or containing folder for processing all fastq files recursively).");
		System.out.println("-o %s \t\tOutput folder");
		System.out.println("-uniqueBarcode\tIf the fastq file(s) contain(s) only one barcode (for e.g. after demultiplexing), this option can be used for searching the specific barcode (most occuring) in the construct and trimming it when present.");
		System.out.println("-polyA %i \tTrim polyA strings that have more than this length (without mismatch), and all 3' string that follows [default=6]");
		System.out.println("-minLength %i \tIf resulting trimmed reads are < this number, it is removed from the output fastq file  [default=10]");
	}
	
	public static void printHelp()
	{
		System.out.println("BRBseq-tools 1.1\n\nOptions:");
		System.out.println("CreateDGEMatrix\tCreate the DGE Matrix (counts + UMI) from R2 aligned BAM and R1 fastq");
		System.out.println("Demultiplex\tCreate demultiplexed fastq files from R1+R2 fastq (use this if you want to process them independently, if not, use the CreateDGEMatrix option)");
		System.out.println("Trim\t\tFor trimming BRBseq construct in R2 fastq file or demultiplexed fastq files.");
		System.out.println("AnnotateBAM\t\tFor annotating the BAM file using UMIs/Barcodes from the fastq file.");
	}
}
