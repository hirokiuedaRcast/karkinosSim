package jp.ac.utokyo.karkinosSim.candidategen;

public class Position {
	
	
	String chr;
	int pos;
	
	public String getChr() {
		return chr;
	}
	public void setChr(String chr) {
		this.chr = chr;
	}
	public int getPos() {
		return pos;
	}
	public void setPos(int pos) {
		this.pos = pos;
	}
	
	boolean indel;
	boolean insersion;

	public boolean isIndel() {
		return indel;
	}
	public void setIndel(boolean indel) {
		this.indel = indel;
	}
	public boolean isInsersion() {
		return insersion;
	}
	public void setInsersion(boolean insersion) {
		this.insersion = insersion;
	}

}
