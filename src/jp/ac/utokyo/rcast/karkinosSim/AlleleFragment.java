package jp.ac.utokyo.rcast.karkinosSim;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class AlleleFragment {
	
	List<SNPInfoBean> snpList = new ArrayList<SNPInfoBean>();
	int start;
	int end;
	int count;
	
	public String getAInfoStr() {
		
		StringBuffer sb = new StringBuffer();
		for(SNPInfoBean info:snpList){
			String s = info.ref+"/"+info.alt+"("+info.modeSnp+")";	
			sb.append(s+"-");
		}
		return sb.toString();
	
	}
	
	public void add(SNPInfoBean bean){
		
		int pos = bean.pos;
		snpList.add(bean);
		if(start==0||pos < start){
			start = pos;
		}
		if(end==0||pos > end){
			end = pos;
		}
		count++;
	}

	Map<Integer, AllelleBase> allelemap;
	public Map<Integer, AllelleBase> getAllelemap() {
		return allelemap;
	}

	public void setAllelemap(Map<Integer, AllelleBase> allelemap) {
		this.allelemap = allelemap;
	}

	public String getALMapStr(){
		
		StringBuffer sb = new StringBuffer();
		Iterator<Integer> ite = allelemap.keySet().iterator();
		while(ite.hasNext()){
			int pos = ite.next();
			AllelleBase ab = allelemap.get(pos);
			sb.append(pos+ ab.alleleA+"|"+ab.alleleB);
		}		
		return sb.toString();
		
	}
	


}
