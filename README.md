![](https://img.shields.io/badge/build-passing-green.svg)
![](https://img.shields.io/badge/version-1.6.1-blue.svg)
![](https://img.shields.io/badge/picard-2.9.0-blue.svg)
![](https://img.shields.io/badge/java-1.8-red.svg)
[![DOI](https://zenodo.org/badge/109837582.svg)](https://zenodo.org/badge/latestdoi/109837582)

# BRB-seq Tools
A suite of tools for the pre-processing of BRB-seq data (bulk RNA-seq)

## :warning: Recent notes
BRB-seq tools suite was created in the early days of multiplexed libraries, when there was not many other alternatives to analyze BRB-seq libraries. Now, this is not the case anymore, so we would recommend using **STARsolo** instead, which should produce similar results in a single command.

In particular, **do NOT use the output of the Trimming step of BRB-seq Tools as input for STARsolo**, as this will not produce correct UMI values (without displaying an error message). STARsolo can provide read trimming that matches BRB-seq specificities using the ``--clipAdapterType CellRanger4`` option

An example of command that you can run to align + demultiplex BRB-seq libraries:
```bash
STAR --runMode alignReads --soloUMIdedup 1MM_Directional --soloCBmatchWLtype 1MM --clipAdapterType CellRanger4 --outSAMmapqUnique 60 --outSAMunmapped Within --soloStrand Forward --quantMode GeneCounts --genomeDir STAR_Index/ --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 16 --soloCellFilter None --soloCBwhitelist lib_example_barcodes.txt --soloFeatures Gene --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/lib_example/ --readFilesIn lib_example_R2.fastq.gz lib_example_R1.fastq.gz
```

> **Note:** Here, the ``lib_example_barcodes.txt`` file should be a list of barcodes, one per line.

Of course, you can also add some other options like:
- ``--runThreadN 12 --outBAMsortingThreadN 12`` to run STAR in parallel on 12 threads (or more, or less)
- ``--soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 16`` this depends on your R1 sequence (check in the fastq file), and the barcodes that you are using

## Download software
BRB-seq command-line tools are provided as a [single executable jar file](../master/releases/BRBseqTools-1.6.1.jar?raw=true).
The .jar file contains all required materials and can be run on any terminal.

## Dependencies
### Java version
For the tools to run properly, you must have Java >=1.8 installed. 

To check your java version, open your terminal application and r
un the following command:

```bash
java -version
```

If the output looks something like java version "1.8.x", you are good to go. 
If not, you may need to update your version; see the [Oracle Java website](http://www.oracle.com/technetwork/java/javase/downloads/) to download the latest JRE (for users) or JDK (for developers).

### Picard
The software relies on [Picard Tools](http://broadinstitute.github.io/picard/), but the Picard JAR is already embedded in the released JAR, so no need to install it yourself.

## Sequencing output
After sequencing your BRB-seq libraries, you should obtain two fastq files per library: 
* R1 fastq file: This file should contain the barcodes and UMIs of the BRBseq construct. Usually, the barcode comes before the UMI sequence. And the UMI sequence is optional.
* R2 fastq file: This file should contain the exact same number of reads (and read names) than the R1 file. Except that the sequence of the reads are the sequences representing the RNA fragments.

BRB-seq tools is a suite dedicated to help you analyze these data, until the generation of the output count/UMI matrix.
For further analyses (filtering, normalization, dimension reduction, clustering, differential expression), we recommend using [ASAP web portal](https://www.ncbi.nlm.nih.gov/pubmed/28541377) that you can freely access at [asap.epfl.ch](https://asap.epfl.ch).

## Aligner
BRB-seqTools was mainly tested with the STAR aligner. So far it seems to not working with the hisat2 aligner. Please contact me if it does not work with any other aligner, I'll try to have a look in the future, if this is a recurrent problem.

## Installation
To check that BRB-seq Tools is working properly, run the following command:

```bash
java -jar BRBseqTools.jar
```
As shown in the output of this command, BRB-seq tool suite allows the user to run 4 possible tools.

```
Demultiplex     Create demultiplexed fastq files from R1+R2 fastq (use this if you want to process them independently, thus doing the next steps without BRBseqTools)
CreateDGEMatrix Create the DGE Matrix (counts + UMI) from R2 aligned BAM and R1 fastq
Trim            For trimming BRBseq construct in R2 fastq file or demultiplexed fastq files
AnnotateBAM     For annotating the R2 BAM file using UMIs/Barcodes from the R1 fastq and/or GTF and/or Barcode files
ExtractReadCountMatrix  For generating a read count matrix from a STAR-aligned BAM file.
```

## Example of full pipeline (using STAR)

First if you don't already have an indexed genome, you need to build the index (for STAR)
```bash
# Download last release of your species of interest .fasta file from Ensembl (or any other database that you'd prefer to use). Here e.g. Homo sapiens from Ensembl
wget ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# Unzip
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# Download last release of homo sapiens .gtf file from Ensembl
wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
# Unzip
gzip -d Homo_sapiens.GRCh38.90.gtf.gz
# Generate the STAR genome index (in 'STAR_Index' folder) => This require ~30G RAM
mkdir STAR_Index/
STAR --runMode genomeGenerate --genomeDir STAR_Index/ --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.90.gtf
# Note: you can add the argument '--runThreadN 4' for running 4 threads in parallel (or more if your computer has enough cores)
# Note: if genome is not human, you need to add/tune the argument '--genomeSAindexNbases xx' with xx = log2(nbBasesInGenome)/2 - 1
#	For e.g. for drosophila melanogaster xx ~ 12.43, thus you should add the argument '--genomeSAindexNbases 12'
```

When the index is built, you will never need to rebuild it. Then you can use only the following script:
```bash
# (Optional) Trim the read containing the sequence fragments (generates a 'lib_example_R2.trimmed.fastq.gz' file)
java -jar BRBseqTools.jar Trim -f lib_example_R2.fastq.gz
# Create output folder
mkdir BAM/
# Align only the R2 fastq file (using STAR, no sorting/indexing is needed)
STAR --runMode alignReads --genomeDir STAR_Index/ --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix BAM/ --readFilesIn lib_example_R2.trimmed.fastq.gz
# Note: you can add the argument '--runThreadN 4' for running 4 threads in parallel (or more if your computer has enough cores)
# Note: the '--outFilterMultimapNmax 1' option is recommended for removing multiple mapping reads from the output BAM
# (optional) Rename the output aligned BAM
mv BAM/Aligned.out.bam BAM/lib_example_R2.bam
# Demultiplex and generate output count/UMI matrix
java -jar BRBseqTools.jar CreateDGEMatrix -f lib_example_R1.fastq.gz -b BAM/lib_example_R2.bam -c lib_example_barcodes.txt -gtf Homo_sapiens.GRCh38.90.gtf -p BU -UMI 14
# Note: This example suppose that R1 has barcode followed by 14bp UMI
# Note: The ‘-p’ option specifies the pattern in the R1 fastq files. If you sequenced only the 6bp barcode then you should use ‘-p B’, but if you sequenced the barcode + 10bp UMIs, then it should be ‘-p BU’ (for Barcode followed by UMI) and then you need to specify the UMI length with ‘-UMI 10’. 
# Note: In case you sequenced some bp you don’t want to use, you can use the ‘?’ character. For e.g. you sequenced 10bp UMI but the R1 contains extra 4bp after the UMI that are not UMI, then you should use “-p BU???? -UMI 10”
# Note: 'lib_example_barcodes.txt' should be created by the user and should contain the mapping between the barcode and the sample name¹
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

Then you can load the generated 'output.dge.reads.txt' count matrix in R and performs your analyses (or 'output.dge.umis.txt' if you prefer working with UMIs).
Or you can upload this file to https://asap.epfl.ch and run the analysis pipeline online.

## Usage
Here follows the description of each tool in more details.

### ExtractReadCountMatrix
This tool was recently added to BRBseq-Tools and is useful when using STARsolo for performing the whole pipeline (instead of BRBseq-Tools). The issue with STARsolo is that it does not generate the read count matrix (only the UMI count matrix). ExtractReadCountMatrix was specifically designed to retrieve the read count matrix from the BAM output from STARsolo.
> **Note:** This tool allows to perform the whole pipeline with one STARsolo command (and eventually trim via cutadapt). It will probably become the method of choice for processing BRB-seq datasets, given STARsolo performances.

Options:
```
-b %s           [Required] Path of STAR-aligned BAM file to analyze.
-c %s           [Required] Path of Barcode/Samplename mapping file.
-gtf %s         [Required] Path of GTF file.
-o %s           Output folder
-chunkSize %i   Maximum number of reads to be stored in RAM (default = 1000000)
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

Example:
```bash
java -jar BRBseqTools.jar ExtractReadCountMatrix -b STAR_solo_aligned.bam -c lib_example_barcodes.txt -gtf hg38.gtf
```

### Demultiplex
This tool is used when you need to generate all the fastq files corresponding to all your multiplexed samples. You end up with one fastq file per sample.
> **Note:** If you use the "AnnotateBAM" tool, you can also demultiplex the annotated BAM file using GATK SplitReads

Options:
```
-r1 %s          [Required] Path of R1 FastQ files (containing barcode and optionally UMIs) [can be gzipped or raw].
-r2 %s          [Required] Path of R2 FastQ files (containing read sequence) [can be gzipped or raw].
-c %s           [Required] Path of Barcode/Samplename mapping file¹.
-n %i           Number of allowed difference with the barcode [Default = 1]. Ambiguous barcodes (same distance from two or more existing barcodes) will be automatically discarded.
-o %s           Output folder
-p %s           Barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).
                        'B' [required] is used for specifying the barcode position.
                        'U' [optional] is used for specifying a UMI value position.
                        Character '?' [optional] can be used to ignore specific nucleotides.
-UMI %i         If your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI.
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

Example:
```bash
java -jar BRBseqTools.jar Demultiplex -r1 lib_example_R1.fastq.gz -r2 lib_example_R2.fastq.gz -c lib_example_barcodes.txt -p BU -UMI 14
```
or, if no UMI:
```bash
java -jar BRBseqTools.jar Demultiplex -r1 lib_example_R1.fastq.gz -r2 lib_example_R2.fastq.gz -c lib_example_barcodes.txt -p B
```

> **Note:** When using this tool, UMIs are kept as indexes in all output .fastq files

### CreateDGEMatrix
This tool is used when you don't need the intermediary .fastq & .bam files from all your multiplexed samples.
It greatly simplifies the demultiplexing and analysis, and directly generates a workable count/UMI matrix.

Options:
```
-f %s           [Required] Path of R1 FastQ file [can be gzipped or raw].
-b %s           [Required] Path of R2 aligned BAM file [do not need to be sorted or indexed].
-c %s           [Required] Path of Barcode/Samplename mapping file¹.
-gtf %s         [Required] Path of GTF file [can be gzipped or raw].
-s %s           Do you want to count only reads falling on same strand than gene? [no, yes, reverse] (default = yes since BRB-seq is stranded protocol)
-n %i           Number of allowed difference with the barcode [ambiguous reads will be automatically discarded].
-o %s           Output folder
-t %s           Path of existing folder for storing temporary files
-chunkSize %i   Maximum number of reads to be stored in RAM (default = 10000000)
-p %s           Barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).
                        'B' [required] is used for specifying the barcode position.
                        'U' [optional] is used for specifying a UMI value position.
                        Character '?' [optional] can be used to ignore specific nucleotides.
-UMI %i         If your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI.
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

Example:
```bash
java -jar BRBseqTools.jar CreateDGEMatrix -f lib_example_R1.fastq.gz -b lib_example_R2.bam -c lib_example_barcodes.txt -gtf Homo_sapiens.GRCh38.90.gtf.gz -p BU -UMI 14
```

> **Note:** The original BRB-seq protocol contains a UMI construct. But you can also generate a library without UMIs, or even not sequence the UMIs from R2 read. If UMIs are present, both UMI and read count matrices will be generated. If not, only the read count table will be generated.

### AnnotateBAM
This tool is used when you need specific downstream analyses of your BAM files that are not handled by BRBseqTools. It creates an annotated BAM file containing the UMI/Barcodes and can also have annotated genes information. 
> **Note:** This tool allows you to use pipelines such as Picard MarkDuplicates, that can use the UMI/Barcode information for better quantifying the duplicated reads (Picard options BARCODE_TAG=BC READ_ONE_BARCODE_TAG=BX), you can also demultiplex the annotated BAM file using GATK SplitReads

Options:
```
-f %s           [Required] Path of R1 FastQ file [can be gzipped or raw].
-b %s           [Required] Path of R2 aligned BAM file [do not need to be sorted or indexed].
-c %s           Path of Barcode/Samplename mapping file¹.
-gtf %s         Path of GTF file [can be gzipped or raw] fir annoting mapped genes
-n %i           Number of allowed difference with the barcode [ambiguous reads will be automatically discarded].
-s %s           Do you want to annotate genes only in case where reads are falling on same strand? [no, yes, reverse] (default = yes since BRB-seq is stranded protocol
-rg             Add Read Group information (for GATK / Picard MarkDuplicates / etc...). Can replace gatk AddOrReplaceReadGroups.
-o %s           Output folder
-t %s           Path of existing folder for storing temporary files
-chunkSize %i   Maximum number of reads to be stored in RAM (default = 10000000)
-p %s           Barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).
                        'B' [required] is used for specifying the barcode position.
                        'U' [optional] is used for specifying a UMI value position.
                        Character '?' [optional] can be used to ignore specific nucleotides.
-UMI %i         If your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI.
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

Example:
```bash
java -jar BRBseqTools.jar AnnotateBAM -f lib_example_R1.fastq.gz -b lib_example_R2.bam -p BU -UMI 14
```

> **Note:** Most of the options are optionals here. If not put, then the corresponding tag will simply not be added to the annotated BAM.

### Trim
This tool is used to trim the BRB-seq construct (SMART oligo + Barcode + UMI) & the polyA sequences from the R2 fastq file.<br />
It can also be used to trim individually all the demultiplexed .fastq files if you used the "Demultiplex" tool (in this case, use the -uniqueBarcode option).<br />
The tool will TRIM the reads (i.e. cut the construct from the read, and write the remaining of the read). Some reads may not be written in case the remaining of the reads is below the minimum allowed length (-minLength option).

Options:
```
-f %s           [Required] Path of FastQ file to trim (or containing folder for processing all fastq files at once).
-o %s           Output folder
-uniqueBarcode  If the fastq file(s) contain(s) only one barcode (for e.g. after demultiplexing), this option can be used for searching the specific barcode (most occuring) in the construct and trimming it when present.
-polyA %i       Trim polyA strings that have more than this length (without mismatch), and all 3' string that follows [default=6]
-minLength %i   If resulting trimmed reads are < this number, it is removed from the output fastq file  [default=10]
```

Example:
```bash
java -jar BRBseqTools.jar Trim -f lib_example_R2.fastq.gz
```

> **Note:** If you use STAR for alignment, this step is optional, as it will not change much the results of the alignment (our tests have shown that the improvement is real but very minor)

## Directory content
* **src**: all source files required for compilation
* **lib**: all JAR dependencies required for compilation / execution
* **releases**: latest release of BRB-seq Tools

## Author
Vincent Gardeux - vincent.gardeux@epfl.ch
