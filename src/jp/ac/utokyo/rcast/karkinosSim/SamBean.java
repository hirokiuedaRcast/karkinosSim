package jp.ac.utokyo.rcast.karkinosSim;

import htsjdk.samtools.SAMRecord;

public class SamBean {

	SAMRecord sam1;
	SAMRecord sam2;

	public SAMRecord getSam1() {
		return sam1;
	}

	public void setSam1(SAMRecord sam1) {
		this.sam1 = sam1;
	}

	public SAMRecord getSam2() {
		return sam2;
	}

	public void setSam2(SAMRecord sam2) {
		this.sam2 = sam2;
	}

	public void add(SAMRecord sam) {

		if (sam.getFirstOfPairFlag()) {
			sam1 = sam;
		} else {
			sam2 = sam;
		}

	}
	
	public int getMin(){
		int smin1 =Integer.MAX_VALUE;
		if(sam1!=null){
			smin1=sam1.getAlignmentStart();
		}
		int smin2 =Integer.MAX_VALUE;
		if(sam2!=null){
			smin2=sam2.getAlignmentStart();
		}
		return Math.min(smin1,smin2);
	}
	public int getMax(){
		int smin1 = 0;
		if(sam1!=null){
			smin1=Math.max(sam1.getAlignmentStart(),sam1.getAlignmentEnd());
		}
		int smin2 = 0;
		if(sam2!=null){
			smin2=Math.max(sam2.getAlignmentStart(),sam2.getAlignmentEnd());
		}
		return Math.max(smin1,smin2);
	}
	
	public boolean bothExsist(){
		return sam1!=null && sam2!=null;
	}

	public boolean contain(int pos) {

		if (containSam(sam1, pos))
			return true;
		if (containSam(sam2, pos))
			return true;
		return false;

	}

	//
	private boolean containSam(SAMRecord sam, int pos) {

		if(sam==null)return false;
		if(sam.getReadUnmappedFlag())return false;
		
		int start =sam.getAlignmentStart();
		int end = sam.getAlignmentEnd();
		
		if((pos>=start)&&(pos<=end)){
			return true;
		}		
		return false;
		
	}

}
