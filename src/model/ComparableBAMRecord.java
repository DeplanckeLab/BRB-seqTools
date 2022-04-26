package model;

import htsjdk.samtools.SAMRecord;

public class ComparableBAMRecord implements Comparable<ComparableBAMRecord>
{
	public SAMRecord record = null;

	public ComparableBAMRecord(SAMRecord record) 
	{
		this.record = record;
	}

	@Override
	public int compareTo(ComparableBAMRecord o) 
	{
		return this.record.getReadName().compareTo(o.record.getReadName());
	}
	
}
