import java.io.File;
import java.util.HashSet;
import java.util.TreeSet;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import model.ComparableBAMRecord;
import model.Parameters;
import model.Read;
import model.TmpFileManager;
import tools.Utils;

public class AnnotateBAMManager 
{
	public static SAMFileHeader header = null;
	
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
		header = samReader.getFileHeader();
		SAMRecordIterator it = samReader.iterator();
		HashSet<String> uniqueGenes = new HashSet<String>();
		// Start reading the BAM file
		TreeSet<ComparableBAMRecord> lines = new TreeSet<ComparableBAMRecord>();
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
			ComparableBAMRecord bamRecord = new ComparableBAMRecord(samRecord);
			Parameters.nbReads++;
			
			if(Parameters.inputGTFFile != null)
			{
				if(samRecord.getReadUnmappedFlag())
				{
					Parameters.unmapped++;
					bamRecord.record.setAttribute("CO", "__not_aligned");
				}
				else if(samRecord.getMappingQuality() < 10)
				{
					Parameters.toolowAqual++;
					bamRecord.record.setAttribute("CO", "__too_low_aQual");
				}
				else
				{
					boolean toProcess = true;
					if(samRecord.getSupplementaryAlignmentFlag())
					{
						Parameters.notUnique++;
						if(!Parameters.keep_multiple_mapped_reads)
						{
							bamRecord.record.setAttribute("CO", "__alignment_not_unique");
							toProcess = false;
						}
					}
					if(toProcess)
					{
						HashSet<String> overlappingGenes = Utils.getOverlappingGenes(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getCigar(), samRecord.getReadNegativeStrandFlag());
						if(overlappingGenes.size() == 0)
						{
							Parameters.noFeature++;
							bamRecord.record.setAttribute("CO", "__no_feature");
						} 
						else if(overlappingGenes.size() == 1)	
						{
							Parameters.mapped++;
							String gene = overlappingGenes.iterator().next();
							uniqueGenes.add(gene);
							bamRecord.record.setAttribute("CO", gene);
						}
						else
						{
							Parameters.ambiguous++;
							bamRecord.record.setAttribute("CO", "__ambiguous");
						}
					}
				}
			}
			lines.add(bamRecord);
			
			if(Parameters.nbReads %Parameters.chunkSize == 0) 
			{
				TmpFileManager.createTmpFileB(header, Parameters.tmpFolder + "sam."+Parameters.nbTmpBAM+".tmp", lines);
				Parameters.nbTmpBAM++;
				lines = new TreeSet<ComparableBAMRecord>();
				System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
		}
		samReader.close();
		
		TmpFileManager.createTmpFileB(header, Parameters.tmpFolder + "sam."+Parameters.nbTmpBAM+".tmp", lines);
		System.out.println(Parameters.nbReads + " reads were processed from BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		System.out.println("Created " + Parameters.nbTmpBAM + " temporary BAM files");
		
		if(Parameters.inputGTFFile != null)
		{
			System.out.println(uniqueGenes.size() + " Unique genes were detected (at least one read).");
			System.out.println(Parameters.mapped + " 'Mapped' reads (" + Parameters.pcFormatter.format(((float)Parameters.mapped / Parameters.nbReads) * 100) + "%)");
			System.out.println(Parameters.ambiguous + " 'Ambiguous' reads (" + Parameters.pcFormatter.format(((float)Parameters.ambiguous / Parameters.nbReads) * 100) + "%)");
			System.out.println(Parameters.noFeature + " 'No Features' reads (" + Parameters.pcFormatter.format(((float)Parameters.noFeature / Parameters.nbReads) * 100) + "%)");
			System.out.println(Parameters.notUnique + " 'Not Unique' reads (" + Parameters.pcFormatter.format(((float)Parameters.notUnique / Parameters.nbReads) * 100) + "%)");
			System.out.println(Parameters.unmapped + " 'Not Aligned' reads (" + Parameters.pcFormatter.format(((float)Parameters.unmapped / Parameters.nbReads) * 100) + "%)");
			System.out.println(Parameters.toolowAqual + " 'Too Low aQual' reads (" + Parameters.pcFormatter.format(((float)Parameters.toolowAqual / Parameters.nbReads) * 100) + "%)");
		}
	}
	
	/**
	 * Read temporary sorted BAM/fastq files
	 * @throws Exception Yes I know...
	 */
	public static void annotate() throws Exception
	{
		Long start = System.currentTimeMillis();
		System.out.println("\nStarting annotation of BAM...");
		TmpFileManager fm_fastq = new TmpFileManager("fastq", Parameters.nbTmpFastq);
		TmpFileManager fm_bam = new TmpFileManager("sam", Parameters.nbTmpBAM);
		
		// Modify the header to add read groups if requested
		if(Parameters.addReadGroup)
		{
			for(String id:Parameters.readGroups)
			{
				SAMReadGroupRecord rgr = new SAMReadGroupRecord(id);
				rgr.setLibrary(Parameters.libname); // Assume that there is only one library in the file
				rgr.setSample(id.split("\\.")[1]);
				rgr.setPlatform(Parameters.platform);
				rgr.setPlatformUnit(id); // Redundant with ID
				header.addReadGroup(rgr);
			}
		}
		
		// Read Tmp files and generate the new BAM on the go
		int nbReadsWritten = 0;
		int lastWritten = nbReadsWritten; // For not writing twice the same line [Mehhhh]
		SAMFileWriterFactory samWriterFactory = new SAMFileWriterFactory();
		SAMFileWriter samWriter = samWriterFactory.makeBAMWriter(header, false, new File(Parameters.outputFolder + Parameters.inputBAMFileR2.getName() + ".annotated.bam"));
		Read r_fq = fm_fastq.getFirstRead();
		while(r_fq != null)
		{
			Read r_bam = fm_bam.getReadS(r_fq.name);
			if(r_bam != null) 
			{
				nbReadsWritten++;
				String sampleName = "UNKNOWN";
				if(Parameters.l1 != -1) 
				{
					r_bam.samRecord.setAttribute("BC", r_fq.barcode);
					r_bam.samRecord.setAttribute("QT", r_fq.qualityBC);
					if(!r_fq.barcode.equals("Unknown")) sampleName = Parameters.mappingBarcodeName.get(r_fq.barcode).replaceAll("\\.", "_"); // Define the name if it exists
				}
				if(Parameters.addReadGroup)
				{
					String[] tokens = r_bam.samRecord.getReadName().split(":"); // TODO what if this does not work? Here the example that works is readName="NB500883:254:HL2K3BGX5:3:13509:21553:18111"
					String id = tokens[2] + "." + sampleName + "." + tokens[3];
					r_bam.samRecord.setAttribute("RG", id); // For each read we add the "@RG:ID field"
				}
				if(Parameters.UMILength != -1)
				{
					r_bam.samRecord.setAttribute("BX", r_fq.UMI); // BX or Maybe should be RX since UMI is not corrected
					r_bam.samRecord.setAttribute("QX", r_fq.qualityUMI);
				}
				samWriter.addAlignment(r_bam.samRecord);
			}
			r_fq = fm_fastq.getFirstRead();
			if(nbReadsWritten %Parameters.chunkSize == 0 && lastWritten != nbReadsWritten) 
			{
				lastWritten = nbReadsWritten;
				System.out.println(nbReadsWritten + " reads were written in annotated BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
		}
		samWriter.close();
	
		// Clean tmp files
		fm_bam.closeSAMReaders();
		TmpFileManager.delete("fastq", Parameters.nbTmpFastq);
		TmpFileManager.delete("sam", Parameters.nbTmpBAM);
		System.out.println(nbReadsWritten + " reads were written in annotated BAM file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
	}
}
