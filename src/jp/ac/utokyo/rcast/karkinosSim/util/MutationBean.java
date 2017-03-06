package jp.ac.utokyo.rcast.karkinosSim.util;

public class MutationBean {

	String chr;
	int pos;
	String ref;
	String alt;
	boolean indel;
	boolean pass;
	boolean pass2;
	String filter;
	String filter2;
	float ratio = 0;

	boolean passcheck;
	boolean pass2check;
	
	boolean exsist;

	public boolean isExsist() {
		return exsist;
	}

	public void setExsist(boolean exsist) {
		this.exsist = exsist;
	}

	public boolean isPasscheck() {
		return passcheck;
	}

	public void setPasscheck(boolean passcheck) {
		this.passcheck = passcheck;
	}

	public boolean isPass2check() {
		return pass2check;
	}

	public void setPass2check(boolean pass2check) {
		this.pass2check = pass2check;
	}

	public MutationBean(String[] stra, int mode) {

		if (mode == -1) {
			chr = stra[0];
			pos = Integer.parseInt(stra[1]);
			
		}else if (mode == 0) {
			chr = stra[0];
			pos = Integer.parseInt(stra[1]);
			ref = stra[3];
			alt = stra[4];
			indel = (ref.length() != alt.length());
			
			boolean snp = (stra[6].trim().equals("snp"));

//			pass = (stra[6].equals("PASS")|| snp);
//			pass2 = stra[8].equals("PASS");
			
			pass = stra[6].equals("PASS");
			pass2 = stra[8].equals("PASS");
			filter = stra[6];
			filter2 = stra[8];
			
			
		} else if(mode==1){
			chr = stra[0];
			pos = Integer.parseInt(stra[1]);
			ref = stra[3];
			alt = stra[4];
			indel = false;
			int n= stra.length;
			pass = stra[n-1].equals("KEEP");
			pass2 = pass;
//			pass2 = stra[8].equals("PASS");
//			filter = stra[6];
		}else if(mode==2){
			chr = stra[0];
			pos = Integer.parseInt(stra[1]);
			ref = stra[2];
			alt = stra[3];
			indel = false;
			int n= stra.length;
			pass = stra[12].equals("Somatic");
		}
	}

	public boolean equals(MutationBean mb) {

		return mb.getChr().equals(chr) && mb.getPos() == pos;

	}

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

	public String getRef() {
		return ref;
	}

	public void setRef(String ref) {
		this.ref = ref;
	}

	public String getAlt() {
		return alt;
	}

	public void setAlt(String alt) {
		this.alt = alt;
	}

	public boolean isIndel() {
		return indel;
	}

	public void setIndel(boolean indel) {
		this.indel = indel;
	}

}
