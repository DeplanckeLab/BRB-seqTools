![](https://img.shields.io/badge/build-passing-green.svg)
![](https://img.shields.io/badge/version-1.0-blue.svg)
![](https://img.shields.io/badge/picard-2.9.0-blue.svg)
![](https://img.shields.io/badge/java-1.8-red.svg)

# BRB-seq Tools 1.0
A suite of tools for the pre-processing of BRB-seq data (bulk RNA-seq)

## Download software
BRB-seq command-line tools are provided as a [single executable jar file](../master/releases/BRBseqTools.1.0.jar?raw=true).
The .jar file contains all required materials and can be run on any terminal.

## Dependencies
### Java version
For the tools to run properly, you must have Java 1.8 installed. 

To check your java version, open your terminal application and run the following command:

```bash
java -version
```

If the output looks something like java version "1.8.x", you are good to go. 
If not, you may need to update your version; see the [Oracle Java website](http://www.oracle.com/technetwork/java/javase/downloads/) to download the latest JRE (for users) or JDK (for developers).

### Picard
The software relies on [Picard Tools](http://broadinstitute.github.io/picard/), but the Picard JAR is already embedded in the released JAR, so no need to install it yourself.

## Usage
After sequencing your BRB-seq libraries, you should obtain two fastq files per library: 
* R1 fastq file: This file should contain the barcodes and UMIs of the BRBseq construct. Usually, the barcode comes before the UMI sequence. And the UMI sequence is optional.
* R2 fastq file: This file should contain the exact same number of reads (and read names) than the R1 file. Except that the sequence of the reads are the sequences representing the RNA fragments.

BRB-seq tools is a suite dedicated to help you analyze these data, until the generation of the output count/UMI matrix.
For further analyses (filtering, normalization, dimension reduction, clustering, differential expression), we recommend using [ASAP web portal](https://www.ncbi.nlm.nih.gov/pubmed/28541377) that you can freely access at [asap.epfl.ch](https://asap.epfl.ch).

You have two options depending on what you aim to do with your RNA-seq data:
* ![](https://img.shields.io/badge/Tool-Demultiplex-blue.svg) Perform demultiplexing before aligning your data to the reference genome. This option will generate one fastq file per sample. Every fastq file can then be aligned and processed independently with standard RNA-seq pipelines.
* ![](https://img.shields.io/badge/Tool-CreateDGEMatrix-blue.svg) Align the R2 file and demultiplex after alignment. This option allow you to perform only one alignment on the whole R2 fastq file. This will generate one BAM file (no need to be sorted) with all your samples in it. Then 'BRB-seq tools' suite can be used to generate the read and UMI count matrices from the R2 BAM and the R1 fastq files.

### Installation
To check that BRB-seq Tools is working properly, run the following command:

```bash
java -jar BRBseqTools.1.0.jar
```
As shown in the output of this command, BRB-seq tool suite allows the user to run 3 possible tools.

```
Demultiplex     Create demultiplexed fastq files from R1+R2 fastq (use this if you want to process them independently, if not, use the CreateDGEMatrix option)
CreateDGEMatrix Create the DGE Matrix (counts + UMI) from R2 aligned BAM and R1 fastq
Trim            For trimming BRBseq construct in R2 fastq file or demultiplexed fastq files.
```

### Demultiplex ![](https://img.shields.io/badge/Tool-Demultiplex-blue.svg)
This tool is used when you need to generate all the fastq files corresponding to all your multiplexed samples. You end up with one fastq file per sample.

Options:
```
-r1 %s          Path of R1 FastQ files (containing barcode and optionally UMIs) [can be gzipped or raw].
-r2 %s          Path of R2 FastQ files (containing read sequence) [can be gzipped or raw].
-c %s           Path of Barcode/Samplename mapping file¹.
-n %i           Number of allowed difference with the barcode [Default = 1]. Ambiguous barcodes (same distance from two or more existing barcodes) will be automatically discarded.
-o %s           Output folder
-p %s           Barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file (default = 'BU' i.e. barcode followed by the UMI).
                        'B' [required] is used for specifying the barcode position.
                        'U' [optional] is used for specifying a UMI value position.
                        Character '?' [optional] can be used to ignore specific nucleotides.
-UMI %i         If your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI.
```

Example:
```bash
java -jar BRBseqTools.1.0.jar Demultiplex -r1 lib_example_R1.fastq.gz -r2 lib_example_R2.fastq.gz -c lib_example_barcodes.txt -p BU -UMI 14
```
or, if no UMI:
```bash
java -jar BRBseqTools.1.0.jar Demultiplex -r1 lib_example_R1.fastq.gz -r2 lib_example_R2.fastq.gz -c lib_example_barcodes.txt -p B
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

> **Note:** When using this tool, UMIs are kept as indexes in all output .fastq files

### CreateDGEMatrix ![](https://img.shields.io/badge/Tool-CreateDGEMatrix-blue.svg)
This tool is used when you don't need the intermediary .fastq & .bam files from all your multiplexed samples.
It greatly simplifies the demultiplexing and analysis, and directly generates a workable count/UMI matrix.

Options:
```
-f %s           Path of R1 FastQ file [can be gzipped or raw].
-b %s           Path of R2 aligned BAM file [do not need to be sorted or indexed].
-c %s           Path of Barcode/Samplename mapping file¹.
-gtf %s         Path of GTF file [can be gzipped or raw].
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

Example:
```bash
java -jar BRBseqTools.1.0.jar CreateDGEMatrix -f lib_example_R1.fastq.gz -b lib_example_R2.bam -c lib_example_barcodes.txt -gtf Homo_sapiens.GRCh38.90.gtf.gz -p BU -UMI 14
```

¹You can download/edit this **[example of barcode/samplename mapping file](../master/examples/lib_example_barcodes.txt)**

> **Note:** The original BRB-seq protocol contains a UMI construct. But UMIs are not yet proven effective for bulk RNA-seq analysis. As such, you can generate a library without UMIs, or even not sequence the UMIs from R2 read. If UMIs are present, both UMI and read count matrices will be generated. If not, only the read count table will be generated.

### Trim ![](https://img.shields.io/badge/Tool-Trim-blue.svg)
This tool is used to trim the BRB-seq construct & the polyA sequences from the R2 fastq file.
It can also be used to trim all the demultiplexed .fastq files after the ![](https://img.shields.io/badge/Tool-Demultiplex-blue.svg) step (in this case, use the -uniqueBarcode option)

Options:
```
-f %s           Path of FastQ file to trim (or containing folder for processing all fastq files at once).
-o %s           Output folder
-uniqueBarcode  If the fastq file(s) contain(s) only one barcode (for e.g. after demultiplexing), this option can be used for searching the specific barcode (most occuring) in the construct and trimming it when present.
-polyA %i       Trim polyA strings that have more than this length (without mismatch), and all 3' string that follows [default=6]
-minLength %i   If resulting trimmed reads are < this number, it is removed from the output fastq file  [default=10]
```

Example:
```bash
java -jar BRBseqTools.1.0.jar Trim -f lib_example_R2.fastq.gz
```

> **Note:** If you use STAR for alignment, this step is optional, as it will not change much the results of the alignment (our tests have shown that the improvement is real but very minor)

## Directory content
* **src**: all source files required for compilation
* **lib**: all JAR dependencies required for compilation / execution
* **releases**: latest release of BRB-seq Tools

## Author
Vincent Gardeux - vincent.gardeux@epfl.ch
