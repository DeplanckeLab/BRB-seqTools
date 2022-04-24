import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import model.Parameters;
import tools.Utils;

public class ExtractReadCountMatrixManager
{
	
	/**
	 * Using Picard to read the reads from the BAM file created by the alignment tool
	 */
	public static void extract()
	{
		Long start = System.currentTimeMillis();
		System.out.println("\nReading the reads from the BAM file...");
		HashSet<String> uniqueGenes = new HashSet<String>();
		int[][] count_matrix = new int[Parameters.barcodeIndex.size()][Parameters.geneIndex.size()];
		
		try
		{
			SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = samReaderFactory.open(Parameters.inputBAMFileR2);
			SAMRecordIterator it = samReader.iterator();
			
			// Start reading the BAM file
			Parameters.nbReads = 0;
			Parameters.unmapped = 0;
			Parameters.notUnique = 0;
			Parameters.ambiguous = 0;
			Parameters.mapped = 0;
			Parameters.noFeature = 0;
			Parameters.toolowAqual = 0;
			Parameters.notDemultiplexed = 0;

			while(it.hasNext())
			{
				SAMRecord samRecord = it.next();
				Parameters.nbReads++;

				String barcode = (String)samRecord.getAttribute("CB");
				if(barcode == null || barcode.equals("-")) // New in STAR, now unknown barcodes are labelled '-' instead of the CB tag being nonexistent.
				{
					barcode = "Unknown";
					Parameters.notDemultiplexed++;
				}
				
				Integer indexBarcode = Parameters.barcodeIndex.get(barcode);
				if(indexBarcode == null)
				{
					System.err.println("ERROR: This barcode " + barcode + " is not in your barcode file. Please check again.");
					System.exit(-1);
				}
				
				if(samRecord.getReadUnmappedFlag())
				{
					Parameters.unmapped++;
					count_matrix[indexBarcode][Parameters.geneIndex.get("__not_aligned")]++;
				}
				else if(samRecord.getMappingQuality() < 10)
				{
					Parameters.toolowAqual++;
					count_matrix[indexBarcode][Parameters.geneIndex.get("__too_low_aQual")]++;
				}
				else
				{
					boolean toProcess = true;
					if(samRecord.getSupplementaryAlignmentFlag())
					{
						Parameters.notUnique++;
						if(!Parameters.keep_multiple_mapped_reads) toProcess = false;
					}
					if(toProcess)
					{
						Parameters.mapped++;
						String gene_id = (String)samRecord.getAttribute("GX");
						String gene_name = (String)samRecord.getAttribute("GN");
						if ((gene_id == null || gene_id == "-") && (gene_name == null || gene_name== "-"))
						{
							Parameters.noFeature++;
							count_matrix[indexBarcode][Parameters.geneIndex.get("__no_feature")]++;
						}
						else if(gene_id != null && gene_id!="-" && gene_name != null && gene_name != "-")
						{
							Integer indexGene = Parameters.geneIndex.get(gene_id);
							if(indexGene == null)
							{
								System.err.println("ERROR: This gene " + gene_id + " is not in your GTF file. Please check again.");
								System.exit(-1);
							}
							count_matrix[indexBarcode][indexGene]++;
						}
						else
						{
							System.err.println("ERROR: Weird?");
							System.exit(-1);
						}
						//Parameters.geneIndex.put("__ambiguous", Parameters.geneIndex.size())
					}
				}
				if(Parameters.nbReads %Parameters.chunkSize == 0) System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			samReader.close();
		}
		catch(IOException ioe)
		{
			System.err.println(ioe.getMessage());
		}

		System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		System.out.println(uniqueGenes.size() + " Unique genes were detected (at least one read).");
		System.out.println(Parameters.mapped + " 'Mapped' reads (" + Parameters.pcFormatter.format(((float)Parameters.mapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.ambiguous + " 'Ambiguous' reads (" + Parameters.pcFormatter.format(((float)Parameters.ambiguous / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.noFeature + " 'No Features' reads (" + Parameters.pcFormatter.format(((float)Parameters.noFeature / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.notUnique + " 'Not Unique' reads (" + Parameters.pcFormatter.format(((float)Parameters.notUnique / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.unmapped + " 'Not Aligned' reads (" + Parameters.pcFormatter.format(((float)Parameters.unmapped / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.toolowAqual + " 'Too Low aQual' reads (" + Parameters.pcFormatter.format(((float)Parameters.toolowAqual / Parameters.nbReads) * 100) + "%)");
		System.out.println(Parameters.notDemultiplexed + " 'Aligned but not demultiplexed' reads (" + Parameters.pcFormatter.format(((float)Parameters.notDemultiplexed / Parameters.nbReads) * 100) + "%)");
		
		// Write final matrices
		createOutputDGE(count_matrix);
	}
	
	/**
	 * Read temporary sorted BAM/fastq files and generate DGE
	 * @throws Exception Yes I know...
	 */
	public static void createOutputDGE(int[][] counts)
	{
		try
		{
			// Create the read count and transcript count matrices
			BufferedWriter bw_reads = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.dge.reads.txt"));
			BufferedWriter bw_reads_detailed = new BufferedWriter(new FileWriter(Parameters.outputFolder + "output.dge.reads.detailed.txt"));
			
			// Sorted indexes
			String[] sortedGeneKeys = Utils.sortKeys(Parameters.geneIndex);
			String[] sortedBarcodeKeys = Utils.sortKeys(Parameters.barcodeIndex);
			
			// Header
			bw_reads.write("Gene_id"); 
			bw_reads_detailed.write("Gene_id\tGene_name"); 
			for(String barcode:sortedBarcodeKeys) 
			{ 
				String mappedBarcode = Parameters.mappingBarcodeName.get(barcode);
				if(mappedBarcode != null) bw_reads.write("\t" + mappedBarcode); 
				if(mappedBarcode == null) mappedBarcode = "Unknown_Barcode";
				bw_reads_detailed.write("\t" + mappedBarcode); 
			}
			bw_reads.write("\n");  
			bw_reads_detailed.write("\n"); 
			
			// Actual values
			for(String gene:sortedGeneKeys)
			{
				HashSet<String> mappedGene = Parameters.mappingGeneIdGeneName.get(gene);
				if(mappedGene != null) mappedGene.remove(gene);
				if(mappedGene != null) bw_reads.write(gene);
				bw_reads_detailed.write(gene + "\t" + Utils.toString(mappedGene)); 
				for(String barcode:sortedBarcodeKeys) 
				{
					String mappedBarcode = Parameters.mappingBarcodeName.get(barcode);
					if(mappedGene != null && mappedBarcode != null) bw_reads.write("\t" + counts[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)]);
					bw_reads_detailed.write("\t" + counts[Parameters.barcodeIndex.get(barcode)][Parameters.geneIndex.get(gene)]);
				}
				if(mappedGene != null) bw_reads.write("\n");
				bw_reads_detailed.write("\n");
			}
			bw_reads.close(); 
			bw_reads_detailed.close();
		}
		catch(IOException ioe)
		{
			System.err.println(ioe.getMessage());
			System.exit(-1);
		}
		System.out.println("DGE count matrix written in " + Parameters.outputFolder);
	}
}
