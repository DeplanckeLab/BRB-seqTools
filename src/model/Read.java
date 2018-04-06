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
	
	public void write(BufferedWriter bw, String UMI) throws IOException
	{
		if(rawData != null && rawData.length == 4)
		{
			if(UMI != null) bw.write(rawData[0] + Parameters.separator + UMI + "\n");
			else bw.write(rawData[0] + "\n");
			bw.write(rawData[1] + "\n");
			bw.write(rawData[2] + "\n");
			bw.write(rawData[3] + "\n");
		}
	}
}
