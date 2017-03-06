package jp.ac.utokyo.rcast.karkinosSim;

public class SNPInfoBean {
	
	public SNPInfoBean(String chr2,int pos2,float af2, String id2, String ref, String alt, int modeSnp) {
		chr = chr2;
		pos = pos2;
		af = af2;
		id = id2;
		
		this.ref = ref; 
		this.alt = alt;
		this.modeSnp = modeSnp;
	}
	String id = null;
	String chr;
	int pos;
	int chainPos;
	float af;	
	String ref;
	String alt;
	int modeSnp;
	
	
	boolean readsSopport;
	
	

	public boolean isIndel() {
		
		boolean noIndel = ref.length()==1&&alt.length()==1;
		return !noIndel;
	}
	
}
