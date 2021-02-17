import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.TreeSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import model.Parameters;
import model.Read;
import model.TmpFileManager;
import model.UMI;
import tools.Utils;

public class DGEMatrixManager 
{
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
		HashSet<String> uniqueGenes = new HashSet<String>();
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
			Parameters.nbReads++;
			
			String toAdd = samRecord.getReadName();
			
			if(samRecord.getReadUnmappedFlag())
			{
				Parameters.unmapped++;
				toAdd += "\t__not_aligned\t";
			}
			else if(samRecord.getMappingQuality() < 10)
			{
				Parameters.toolowAqual++;
				toAdd += "\t__too_low_aQual\t";
			}
			else
			{
				boolean toProcess = true;
				if(samRecord.getSupplementaryAlignmentFlag())
				{
					Parameters.notUnique++;
					if(!Parameters.keep_multiple_mapped_reads)
					{
						toAdd += "\t__alignment_not_unique\t";
						toProcess = false;
					}
				}
				if(toProcess)
				{
					HashSet<String> overlappingGenes = Utils.getOverlappingGenes(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag());
					if(overlappingGenes.size() == 0)
					{
						Parameters.noFeature++;
						toAdd += "\t__no_feature\t";
					} 
					else if(overlappingGenes.size() == 1)	
					{
						Parameters.mapped++;
						String gene = overlappingGenes.iterator().next();
						uniqueGenes.add(gene);
						toAdd += "\t" + gene+ "\t";
					}
					else
					{
						Parameters.ambiguous++;
						toAdd += "\t__ambiguous\t";
					}
				}
			}
			lines.add(toAdd + samRecord.getAlignmentStart() + "\t" + samRecord.getAlignmentEnd());

			if(Parameters.nbReads %Parameters.chunkSize == 0) 
			{
				TmpFileManager.createTmpFile(Parameters.tmpFolder + "bam."+Parameters.nbTmpBAM+".tmp", lines);
				Parameters.nbTmpBAM++;
				lines = new TreeSet<String>();
				System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
		}
		samReader.close();
			
		TmpFileManager.createTmpFile(Parameters.tmpFolder + "bam."+Parameters.nbTmpBAM+".tmp", lines);
		System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		System.out.println("Created " + Parameters.nbTmpBAM + " temporary BAM files");
		
		System.out.println(uniqueGenes.size() + " Unique genes were detected (at least one read).");
		System.out.println(Parameters.mapped + " 'Mapped' reads (" + Parameters.pcFormatter.format(((float)Parameters.mapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.ambiguous + " 'Ambiguous' reads (" + Parameters.pcFormatter.format(((float)Parameters.ambiguous / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.noFeature + " 'No Features' reads (" + Parameters.pcFormatter.format(((float)Parameters.noFeature / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.notUnique + " 'Not Unique' reads (" + Parameters.pcFormatter.format(((float)Parameters.notUnique / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.unmapped + " 'Not Aligned' reads (" + Parameters.pcFormatter.format(((float)Parameters.unmapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.toolowAqual + " 'Too Low aQual' reads (" + Parameters.pcFormatter.format(((float)Parameters.toolowAqual / Parameters.nbReads) * 100) + "%)");
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
		
		int nbReadsWritten = 0;
		int lastWritten = nbReadsWritten; // For not writing twice the same line [Mehhhh]
		
		Read r_fq = fm_fastq.getFirstRead();
		while(r_fq != null)
		{
			Read r_bam = fm_bam.getReadB(r_fq.name);
			if(r_bam != null) 
			{
				nbReadsWritten++;
				int x = Parameters.barcodeIndex.get(r_fq.barcode);
				int y = Parameters.geneIndex.get(r_bam.gene);
				counts[x][y]++;
				umis[x][y].addUMI(r_fq.UMI + ":" + r_bam.start + ":" + r_bam.end);
			}
			r_fq = fm_fastq.getFirstRead();
			if(nbReadsWritten %Parameters.chunkSize == 0 && lastWritten != nbReadsWritten) 
			{
				lastWritten = nbReadsWritten;
				System.out.println(nbReadsWritten + " reads were demultiplexed [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
		}

		// Clean tmp files
		TmpFileManager.delete("fastq", Parameters.nbTmpFastq);
		TmpFileManager.delete("bam", Parameters.nbTmpBAM);
		System.out.println(nbReadsWritten + " reads were demultiplexed [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]\n");
		
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
		
		System.out.println("\nStarting writing output file (" + ((Parameters.hammingDistanceUMI == 0)?"without sequencing error correction":"with sequencing error correction: hamming distance <= " + Parameters.hammingDistanceUMI) + ")...");
		String[] sortedKeys = Utils.sortKeys(Parameters.geneIndex);
		int nbUMIs = 0;
		for(String gene:sortedKeys)
		{
			HashSet<String> mappedGene = Parameters.mappingGeneIdGeneName.get(gene);
			if(mappedGene != null) mappedGene.remove(gene);
			if(mappedGene != null) {bw_reads.write(gene); if(Parameters.UMILength != -1) bw_umis.write(gene); }
			bw_reads_detailed.write(gene + "\t" + Utils.toString(mappedGene)); 
			if(Parameters.UMILength != -1) bw_umis_detailed.write(gene + "\t" + Utils.toString(mappedGene));
			for(String barcode:Parameters.barcodeIndex.keySet()) 
			{
				String mappedBarcode = Parameters.mappingBarcodeName.get(barcode);
				nbUMIs += umis[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)].getCorrectedSize();
				if(mappedGene != null && mappedBarcode != null)
				{
					bw_reads.write("\t" + counts[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)]);
					if(Parameters.UMILength != -1) bw_umis.write("\t" + umis[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)].getCorrectedSize());
				}
				bw_reads_detailed.write("\t" + counts[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)]);
				if(Parameters.UMILength != -1) bw_umis_detailed.write("\t" + umis[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)].getCorrectedSize());
			}
			if(mappedGene != null) { bw_reads.write("\n"); if(Parameters.UMILength != -1) bw_umis.write("\n"); }
			bw_reads_detailed.write("\n"); if(Parameters.UMILength != -1) bw_umis_detailed.write("\n");
		}
		System.out.println(nbUMIs + " UMI counts were written (" + (nbReadsWritten - nbUMIs) + " duplicates = "+ Parameters.pcFormatter.format(((nbReadsWritten - nbUMIs) / (float)nbReadsWritten)*100) + "%)");
		bw_reads.close(); bw_reads_detailed.close();
		if(Parameters.UMILength != -1) { bw_umis.close(); bw_umis_detailed.close();}
		System.out.println("DGE count & UMI matrices were created [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
	}
}
