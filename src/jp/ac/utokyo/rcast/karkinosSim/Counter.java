package jp.ac.utokyo.rcast.karkinosSim;

public class Counter {
	
	int n = 1;
	public Counter(String key) {
		this.key = key;
	}
	public Counter() {
		// TODO Auto-generated constructor stub
	}
	public String getKey() {
		return key;
	}
	String key;
	public void inc(){
		n++;
	}
	public int getN(){
		return n;
	}
}
