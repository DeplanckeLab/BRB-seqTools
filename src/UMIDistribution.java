import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import model.Parameters;
import tools.Utils;

public class UMIDistribution 
{
	public static boolean flag = false;
	public static int[][] distribution = null; // [pos][A/C/G/T]
	
	public static void main(String[] args) throws Exception
	{
		Parameters.inputFastQFileR1 = new File("C:/Users/gardeux/Desktop/D1T0A_R1.fastq.gz");
		Parameters.UMIRange[0] = 6;
		Parameters.UMIRange[1] = 16;
		Parameters.outputFolder = "C:/Users/gardeux/Dropbox/SCRBseq/Manuscript/Figures/FigureX/UMI/Distribution/";
		distribution = new int[Parameters.UMIRange[1] - Parameters.UMIRange[0]][20];
		readFastq();
	}
	
	/**
	 *  Reading reads barcodes/UMI from the R1 fastq file to map the UMI/barcode with the read name (lost after alignment)
	 * @throws Exception Yes I know...
	 */
	public static void readFastq() throws Exception
	{
		System.out.println("Reading UMIs from the fastq file...");
		BufferedReader br = Utils.readFastq(Parameters.inputFastQFileR1);
		Set<UMICount> umis = new HashSet<UMICount>();
		Long start = System.currentTimeMillis();
		String umi = nextUMI(br);
		while(umi != null)
		{
			Parameters.nbReads++;
			int i = 0;
			for(char c:umi.toCharArray())
			{
				distribution[i][c - 65]++; // 65 = A, 67 = C, 71 = G, 84 = T, 78 = N
				i++;
			}
			if(Parameters.nbReads %Parameters.chunkSize == 0) 
			{
				System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
			}
			flag = false;
			if(!umis.add(new UMICount(umi, 1)))
			{
				if(!flag) throw new Exception("It failed");
			}
			umi = nextUMI(br);
		}
		br.close();
		
		System.out.println(Parameters.nbReads + " reads were processed from fastq file [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
		System.out.println(umis.size() + " unique UMIs were found");
			
		System.out.println("\nSorting UMIs and generating output file.");
		start = System.currentTimeMillis();
		ArrayList<UMICount> list = new ArrayList<>(umis);
		Collections.sort(list);
		System.out.println("UMIs are sorted [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");

		int sanityCheck = 0;
		BufferedWriter bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.inputFastQFileR1.getName().substring(0, Parameters.inputFastQFileR1.getName().lastIndexOf(".fastq.gz")) + ".unique.umis.count.txt"));
		for(UMICount u:list) 
		{
			sanityCheck += u.count;
			bw.write(u.uniqueUMI + "\t" + u.count + "\t" + Parameters.pcFormatter.format((u.count / (float)Parameters.nbReads) * 100) + "\n");
		}
		bw.close();
		System.out.println("SanityCheck: " + sanityCheck);
		
		bw = new BufferedWriter(new FileWriter(Parameters.outputFolder + Parameters.inputFastQFileR1.getName().substring(0, Parameters.inputFastQFileR1.getName().lastIndexOf(".fastq.gz")) + ".umis.distribution.txt"));
		for(int i = 0; i < distribution.length; i++) bw.write("\t" + (i+1));
		bw.write("\n");
		int[] bases = new int[] {0, 2, 6, 19};//65 = A, 67 = C, 71 = G, 84 = T
		for(int base:bases)
		{
			bw.write((char)(base + 65));
			for(int i = 0; i < distribution.length; i++)
			{
				bw.write("\t" + distribution[i][base]);
			}
			bw.write("\n");
		}
		bw.close();
		
		System.out.println("UMI count distribution is generated [" + Utils.toReadableTime(System.currentTimeMillis() - start) + "]");
	}
	
	/**
	 * Read FastQ file containing the mapping between read name, barcode and UMI.
	 * @param fastQ input file
	 * @return BufferedReader handle
	 * @throws Exception
	 */
	public static String nextUMI(BufferedReader br) throws Exception
	{
		String res = null;
		String header = br.readLine(); // Get first line = READNAME + INDEX
		if(header == null) return null;
		if(!header.startsWith("@")) throw new Exception("Fastq file has formatting issues");
		String indexes = br.readLine().trim(); // Get second line = READ
		try
		{
			res = indexes.substring(Parameters.UMIRange[0], Parameters.UMIRange[1]);
		}
		catch(IndexOutOfBoundsException boe)
		{
			br.close();
			System.err.println("Error while parsing FastQ file");
			System.err.println(boe.getMessage());
			System.exit(-1);
		}
		br.readLine(); // Third line = we don't care
		br.readLine(); // Fourth line = we don't care
		return res;
	}
}

class UMICount implements Comparable<UMICount>
{
	String uniqueUMI = null;
	int count = 0;
	
	public UMICount(String s, int c) 
	{
		uniqueUMI = s;
		count = c;
	}
	
	@Override
	public boolean equals(Object o) // Ugly way of corrupting the array :D
	{
		if(!(o instanceof UMICount)) return false;
		UMICount umi = (UMICount)o;
		if(this.uniqueUMI.equals(umi.uniqueUMI))
		{
			UMIDistribution.flag = true;
			this.count += umi.count;
			umi.count = this.count;
			return true;
		}
		return false;
	}
	
	@Override
	public int hashCode() {
	    return uniqueUMI.hashCode();
	}
	
	@Override
	public int compareTo(UMICount o) 
	{
		if(this.count < o.count) return 1;
		if(this.count > o.count) return -1;
		return 0;
	}
}