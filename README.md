![](https://img.shields.io/badge/build-passing-green.svg)
![](https://img.shields.io/badge/version-1.0-blue.svg)
![](https://img.shields.io/badge/picard-2.9.0-blue.svg)
![](https://img.shields.io/badge/java-1.8-red.svg)

# BRB-seq Tools
A suite of tools for the pre-processing of BRB-seq data (bulk RNA-seq)

## Download software

BRB-seq command-line tools are provided as a single executable jar file. 
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
The software relies on [Picard Tools](http://broadinstitute.github.io/picard/), but the Picard jar is embedded in the released JAR, so no need to install it yourself.

### Usage

To check that BRB-seq Tools is working properly, run the following command:

```bash
java -jar BRBseqTools.jar
```