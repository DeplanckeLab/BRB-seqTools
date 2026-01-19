package model;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

public class Read implements Comparable<Read>
{
	public String name;
	public String barcode;
	public String qualityBC;
	public String UMI;
	public String qualityUMI;
	public String gene;
	public SAMRecord samRecord;
	public String start;
	public String end;
	public String[] rawData; // 4 lines
	public boolean barcodeMatch = false;
	
	@Override
	public int compareTo(Read r2) 
	{
		return this.name.compareTo(r2.name);
	}
	
	public void write(BufferedWriter bw, String barcode, String UMI) throws IOException
	{
		if(rawData != null && rawData.length == 4)
		{
			if(UMI != null)
			{
				bw.write(rawData[0].replaceFirst(name, name + "_" + barcode + "_" + UMI)); // Barcode should be corrected, if possible
			}
			else bw.write(rawData[0]); // In this case I don't even write the barcode
			bw.newLine();
			bw.write(rawData[1]);
			bw.newLine();
			bw.write(rawData[2]);
			bw.newLine();
			bw.write(rawData[3]);
			bw.newLine();
		}
	}
}
