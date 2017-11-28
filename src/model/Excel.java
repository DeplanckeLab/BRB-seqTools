package model;

import java.io.BufferedWriter;
import java.io.IOException;

public class Excel
{
	public int nbReads = 0;
	public int nbRemaininingReads = 0;
	public int nbContaminated = 0;
	public int nbPolyATrimmed = 0;
	public int nbRemoved = 0;
	public int nbPolyABefore = 0;
	
	public void write(BufferedWriter bw) throws IOException
	{
		bw.write("\t"+nbReads+"\t"+nbRemaininingReads+"\t"+nbContaminated+"\t"+nbPolyATrimmed+"\t"+nbRemoved+"\t"+nbPolyABefore+"\n");
	}
	
	public String toString()
	{
		String out = "";
		out += nbReads + " reads in file.\n";
		out += nbRemaininingReads + " are written in the trimmed file.\n";
		out += nbContaminated + " reads were contaminated.\n";
		out += nbRemoved + " reads were removed.\n";
		out += nbPolyATrimmed + " reads were trimmed for polyA (>" + Parameters.polyALength + " 'A's).\n";
		out += "Nb reads with polyA >= " + Parameters.polyALength + " 'A's BEFORE trimming = " + nbPolyABefore;
		return out;
	}
}
