import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import model.GTF;
import model.Parameters;
import model.Read;
import model.TmpFileManager;
import model.UMI;
import tools.Utils;

public class DGEMatrixManager 
{
	/**
	 *  Reading reads barcodes/UMI from the R1 fastq file to map the UMI/barcode with the read name (lost after alignment)
	 * @throws Exception Yes I know...
	 */
	public static void readR1Fastq() throws Exception
	{
		System.out.println("\nReading reads barcodes/UMI from the R1 fastq file...");
		BufferedReader br = Utils.readFastq(Parameters.inputFastQFileR1);
		TreeSet<String> lines = new TreeSet<String>();
		int noBarcodeMatch = 0;
		Long start = System.currentTimeMillis();
		Parameters.nbTmpFastq = 1; // Number of temporary files
		Read read = Utils.nextRead(br);
		while(read != null)
		{
			Parameters.nbReads++;
			if(Parameters.nbReads %Parameters.chunkSize == 0) 
			{
				BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.tmpFolder + "fastq."+Parameters.nbTmpFastq+".tmp"));
				for(String l:lines) bw.write(l + "\n");
				bw.close();
				Parameters.nbTmpFastq++;
				lines = new TreeSet<String>();
				System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			if(!read.barcodeMatch) { noBarcodeMatch++; read.barcode = "Unknown";}
			lines.add(read.name + "\t" + read.barcode + "\t" + read.UMI);
			read = Utils.nextRead(br);
		}
		br.close();
		
		System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.tmpFolder + "fastq."+Parameters.nbTmpFastq+".tmp"));
		for(String l:lines) bw.write(l + "\n");
		bw.close();
		
		System.out.println("Created " + Parameters.nbTmpFastq + " temporary fastq files");
		
		System.out.println(noBarcodeMatch + " reads have no matching barcodes (" + Parameters.pcFormatter.format(((float)noBarcodeMatch / Parameters.nbReads) * 100) + "%)");
	}
	
	/**
	 * Using Picard to read the reads from the BAM file created by the alignment tool
	 * @throws Exception Yes I know...
	 */
	public static void readR2BAM() throws Exception
	{
		Long start = System.currentTimeMillis();
		System.out.println("\nReading the reads from the BAM file...");
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader samReader = samReaderFactory.open(Parameters.inputBAMFileR2);
		SAMRecordIterator it = samReader.iterator();
		// Start reading the BAM file
		TreeSet<String> lines = new TreeSet<String>();
		Parameters.nbReads = 0;
		Parameters.unmapped = 0;
		Parameters.notUnique = 0;
		Parameters.ambiguous = 0;
		Parameters.mapped = 0;
		Parameters.noFeature = 0;
		Parameters.toolowAqual = 0;
		Parameters.nbTmpBAM = 1;
		while(it.hasNext())
		{
			SAMRecord samRecord = it.next();
			if(samRecord.getSupplementaryAlignmentFlag())
			{
				Parameters.notUnique++;
				lines.add(samRecord.getReadName() + "\t__alignment_not_unique");
			}
			else if(!samRecord.getReadUnmappedFlag())
			{
				Parameters.nbReads++;
				HashSet<String> overlappingGenes = getOverlappingGenes(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar());
				if(overlappingGenes.size() == 0)
				{
					Parameters.noFeature++;
					lines.add(samRecord.getReadName() + "\t__no_feature");
				} 
				else if(overlappingGenes.size() == 1)	
				{
					Parameters.mapped++;
					lines.add(samRecord.getReadName() + "\t" + overlappingGenes.iterator().next());
				}
				else
				{
					Parameters.ambiguous++;
					lines.add(samRecord.getReadName() + "\t__ambiguous");
				}
			} 
			else 
			{
				Parameters.unmapped++;
				lines.add(samRecord.getReadName() + "\t__not_aligned");
			}
			if(Parameters.nbReads %Parameters.chunkSize == 0) 
			{
				BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.tmpFolder + "bam."+Parameters.nbTmpBAM+".tmp"));
				for(String l:lines) bw.write(l + "\n");
				bw.close();
				Parameters.nbTmpBAM++;
				lines = new TreeSet<String>();
				System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
		}
		samReader.close();
		
		System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.tmpFolder + "bam."+Parameters.nbTmpBAM+".tmp"));
		for(String l:lines) bw.write(l + "\n");
		bw.close();
		System.out.println("Created " + Parameters.nbTmpBAM + " temporary BAM files");
		
		System.out.println(Parameters.mapped + " 'Mapped' reads (" + Parameters.pcFormatter.format(((float)Parameters.mapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.ambiguous + " 'Ambiguous' reads (" + Parameters.pcFormatter.format(((float)Parameters.ambiguous / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.noFeature + " 'No Features' reads (" + Parameters.pcFormatter.format(((float)Parameters.noFeature / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.notUnique + " 'Not Unique' reads (" + Parameters.pcFormatter.format(((float)Parameters.notUnique / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.unmapped + " 'Not Aligned' reads (" + Parameters.pcFormatter.format(((float)Parameters.unmapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.toolowAqual + " 'Too Low aQual' reads (" + Parameters.pcFormatter.format(((float)Parameters.toolowAqual / Parameters.nbReads) * 100) + "%)");
	}
	
	private static HashSet<String> getOverlappingGenes(String chr, int start, int end, Cigar c) throws Exception
	{
		HashSet<String> res = new HashSet<>();
		List<CigarElement> l = c.getCigarElements();
		int s = start;
		for(CigarElement cigar:l)
		{
			switch(cigar.getOperator())
			{
				case M:
					res.addAll(GTF.findOverlappingGenes(chr, s, s + cigar.getLength()));
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
	 * Read temporary sorted BAM/fastq files and generate DGE
	 * @throws Exception Yes I know...
	 */
	public static void createOutputDGE() throws Exception
	{
		Long start = System.currentTimeMillis();
		System.out.println("\nStarting demultiplexing...");
		TmpFileManager fm_fastq = new TmpFileManager("fastq", Parameters.nbTmpFastq);
		TmpFileManager fm_bam = new TmpFileManager("bam", Parameters.nbTmpBAM);
		
		int[][] counts = new int[Parameters.barcodeIndex.size()][Parameters.geneIndex.size()]; // Sample / gene
		UMI[][] umis = new UMI[Parameters.barcodeIndex.size()][Parameters.geneIndex.size()]; // Sample / gene
		for(int i = 0; i < umis.length; i++) for(int j = 0; j < umis[i].length; j++) umis[i][j] = new UMI();
		
		Read r_fq = fm_fastq.getFirstRead();
		while(r_fq != null)
		{
			Read r_bam = fm_bam.getRead(r_fq.name);
			if(r_bam != null) 
			{
				int x = Parameters.barcodeIndex.get(r_fq.barcode);
				int y = Parameters.geneIndex.get(r_bam.gene);
				counts[x][y]++;
				umis[x][y].umis.add(r_fq.UMI);
			}
			r_fq = fm_fastq.getFirstRead();
		}

		// Clean tmp files
		TmpFileManager.delete("fastq", Parameters.nbTmpFastq);
		TmpFileManager.delete("bam", Parameters.nbTmpBAM);
		
		// Create the read count and transcript count matrices
		BufferedWriter bw_reads = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.dge.reads.txt"));
		BufferedWriter bw_umis = null;
		if(Parameters.UMILength != -1) bw_umis = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.dge.umis.txt"));
		BufferedWriter bw_reads_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.dge.reads.detailed.txt"));
		BufferedWriter bw_umis_detailed = null;
		if(Parameters.UMILength != -1) bw_umis_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.dge.umis.detailed.txt"));
		bw_reads.write("Gene_id"); 
		if(Parameters.UMILength != -1) bw_umis.write("Gene_id");
		bw_reads_detailed.write("Gene_id\tGene_name"); 
		if(Parameters.UMILength != -1) bw_umis_detailed.write("Gene_id\tGene_name");
		for(String barcode:Parameters.barcodeIndex.keySet()) 
		{ 
			String mappedBarcode = Parameters.mappingBarcodeName.get(barcode);
			if(mappedBarcode != null) bw_reads.write("\t" + mappedBarcode); 
			if(mappedBarcode != null && Parameters.UMILength != -1) bw_umis.write("\t" + mappedBarcode);
			if(mappedBarcode == null) mappedBarcode = "Unknown_Barcode";
			bw_reads_detailed.write("\t" + mappedBarcode); 
			if(Parameters.UMILength != -1) bw_umis_detailed.write("\t" + mappedBarcode); 
		}
		bw_reads.write("\n");  bw_reads_detailed.write("\n"); 
		if(Parameters.UMILength != -1) 
		{
			bw_umis.write("\n");
			bw_umis_detailed.write("\n");
		}
		
		String[] sortedKeys = Utils.sortKeys(Parameters.geneIndex);
		for(String gene:sortedKeys)
		{
			String mappedGene = Parameters.mappingGeneIdGeneName.get(gene);
			if(mappedGene == null) mappedGene = "";
			if(!mappedGene.equals("")) {bw_reads.write(gene); if(Parameters.UMILength != -1) bw_umis.write(gene); }
			bw_reads_detailed.write(gene + "\t" + mappedGene); 
			if(Parameters.UMILength != -1) bw_umis_detailed.write(gene + "\t" + mappedGene);
			for(String barcode:Parameters.barcodeIndex.keySet()) 
			{
				String mappedBarcode = Parameters.mappingBarcodeName.get(barcode);
				if(!mappedGene.equals("") && mappedBarcode != null)
				{
					bw_reads.write("\t" + counts[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)]);
					if(Parameters.UMILength != -1) bw_umis.write("\t" + umis[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)].umis.size());
				}
				bw_reads_detailed.write("\t" + counts[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)]);
				if(Parameters.UMILength != -1) bw_umis_detailed.write("\t" + umis[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)].umis.size());
			}
			if(!mappedGene.equals("")) { bw_reads.write("\n"); if(Parameters.UMILength != -1) bw_umis.write("\n"); }
			bw_reads_detailed.write("\n"); if(Parameters.UMILength != -1) bw_umis_detailed.write("\n");
		}
		bw_reads.close(); bw_reads_detailed.close();
		if(Parameters.UMILength != -1) { bw_umis.close(); bw_umis_detailed.close();}
		System.out.println("DGE count & UMI matrices were created [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
	}
}
