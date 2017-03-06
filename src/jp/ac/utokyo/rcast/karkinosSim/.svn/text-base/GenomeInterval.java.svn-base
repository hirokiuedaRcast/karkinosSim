package jp.ac.utokyo.rcast.karkinosSim;

public class GenomeInterval {
	
	String chr;
	int start;
	int end;

	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public boolean include(String chr2, int pos) {
		
		if(chr2.equals(chr)){
			if(pos>=start&&pos<=end)return true;
		}
		return false;
	}
	
	

}
