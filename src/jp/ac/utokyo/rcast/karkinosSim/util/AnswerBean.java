package jp.ac.utokyo.rcast.karkinosSim.util;

import java.util.Map;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinosSim.Counter;

public class AnswerBean {

	Map<String,Counter> amc;
	public Map<String, Counter> getAmc() {
		return amc;
	}
	public void setAmc(Map<String, Counter> amc) {
		this.amc = amc;
	}
	public Map<String, Counter> getMc() {
		return mc;
	}
	public void setMc(Map<String, Counter> mc) {
		this.mc = mc;
	}
	Map<String,Counter> mc;
	
}
