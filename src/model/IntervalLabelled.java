package model;

import htsjdk.tribble.index.interval.Interval;

public class IntervalLabelled extends Interval
{
	public String gene;
	public boolean readNegativeStrandFlag;
	
	public IntervalLabelled(int start, int end, String gene, boolean readNegativeStrandFlag) 
	{
		super(start, end);
		this.gene = gene;
		this.readNegativeStrandFlag = readNegativeStrandFlag;
	}
	
	@Override
	public String toString() {
		return super.toString() + " : " + gene;
	}
}
