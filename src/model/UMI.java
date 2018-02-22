package model;

import java.util.HashSet;

public class UMI 
{
	public HashSet<String> umis = null;
	
	public UMI() 
	{
		umis = new HashSet<>();
	}
	
	public boolean addUMI(String umi) // TODO Should check for sequencing errors mismatches
	{
		return umis.add(umi); // return false if duplicated
	}
}
