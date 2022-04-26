package model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.util.StringUtil;

public class UMI 
{
	public HashSet<String> umis = null;
	public int correctedSize = -1;
	
	public UMI() 
	{
		umis = new HashSet<String>();
	}
	
	public void addUMI(String umi)
	{
		umis.add(umi);
	}
	
	public int getCorrectedSize() throws Exception
	{
		if(correctedSize == -1)
		{
			if(Parameters.hammingDistanceUMI == 0) { correctedSize = umis.size(); return correctedSize; } // If not sequencing error correction
			correctedSize = 0;
			HashMap<String, ArrayList<String>> positions = new HashMap<String, ArrayList<String>>();
			for(String umi:umis)
			{
				int pos = umi.indexOf(":");
				if(pos == -1) throw new Exception("This should not happen");
				String position = umi.substring(pos + 1, umi.length()); // Extract position from String
				ArrayList<String> list = positions.get(position);
				if(list == null) list = new ArrayList<String>();
				list.add(umi.substring(0, pos)); // Extract UMI from String / Should all be different
				positions.put(position, list);
			}
			
			for(String pos:positions.keySet())
			{
				ArrayList<String> list = positions.get(pos);
				if(list.size() > 1) // UMI correction needed
				{
				    for (int i = 0; i < list.size(); i++) 
				    {
				    	boolean flag = false;
				    	for (int j = i + 1; j < list.size(); j++)
				    	{
				    		if (StringUtil.isWithinHammingDistance(list.get(i), list.get(j), Parameters.hammingDistanceUMI))
				    		{
				    			flag = true;
				    			break;
				    		}
				    	}
				    	if(!flag) correctedSize++;
			        }
				}
				else correctedSize += 1;
			}
		}
		return correctedSize; // Should be corrected for sequencing errors mismatches*/
	}
}
