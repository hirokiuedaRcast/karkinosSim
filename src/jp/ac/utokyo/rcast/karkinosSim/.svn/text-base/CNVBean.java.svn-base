package jp.ac.utokyo.rcast.karkinosSim;

public class CNVBean extends GenomeInterval{
	
	
	//cn0-4,7
	int cn; 
	int cnA;
	int cnB;
	boolean focal;
	public CNVBean(String line) {
		
		String[] data = line.split("\t");
		chr = data[0];
		start = toInt(data[1]);
		end = toInt(data[2]);
		cn = toInt(data[3]);
		cnA = toInt(data[4]);
		cnB = toInt(data[5]);
		
		if(data.length>6){
			focal = data[6].equalsIgnoreCase("true");
		}
	}
	public CNVBean() {
		cn = 2;
		cnA = 1;
		cnB = 1;
	}
	private int toInt(String s) {
		
		try{
			return Integer.parseInt(s);
		}catch(Exception ex){
			
		}
		return 0;
	}
	public int getCn() {
		return cn;
	}
	public void setCn(int cn) {
		this.cn = cn;
	}	

}
