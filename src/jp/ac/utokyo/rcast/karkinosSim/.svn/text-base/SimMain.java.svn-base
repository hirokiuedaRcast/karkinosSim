package jp.ac.utokyo.rcast.karkinosSim;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import jp.ac.utokyo.rcast.karkinosSim.bean.SimBean;

public class SimMain {

	public static void main(String[] arg) throws IOException {

		int buffersize = 1000000;

		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		String targetRegion = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		
		Map<String, Integer> map = getChrPortion(buffersize,twobitref);

		Iterator<String> ite = map.keySet().iterator();
		int total = 0;
		while (ite.hasNext()) {

			String key = ite.next();
			int len = map.get(key);
			System.out.println(key+"\t"+len);
			total = total + len;
		}
				
		// get SNP list
		List<SimBean> simList = new ArrayList<SimBean>();

		// while(simList.size()<=buffersize){
		//
		// //SimBean = genPos();
		//
		// }
		//
		// exec();

	}
	
	public static Map<String, Integer> getChrPortion(int buffersize,String twobitref) throws IOException{
		
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		Map<String, Integer> readSize = tgr.getReadSizes(new File(twobitref));
		Iterator<String> ite = readSize.keySet().iterator();
		
		long total = 0;

		while (ite.hasNext()) {

			String key = ite.next();
			int length = readSize.get(key);
			if (length > 10000000) {				
				total = total + length;
			}
		}

		ite = readSize.keySet().iterator();
		Map<String, Integer> readSize2 = new LinkedHashMap<String, Integer>();
		int total2 = 0;
		while (ite.hasNext()) {

			String key = ite.next();
			int length = readSize.get(key);
			if (length > 10000000) {
				int frac = (int)(buffersize * (double)length/(double)total);
				readSize2.put(key, frac);
				total2 = total2+frac;
			}
		}
		if(total2!=buffersize){
			
			int diff = buffersize-total2;	
			String firstkey = readSize2.keySet().iterator().next();
			int val = readSize2.get(firstkey)+diff;
			readSize2.put(firstkey, val);
		}
		return readSize2;
		
	}
	

	public static void exec() {

		//
		double mean = 160;
		double sd = 60;
		int seqlen = 100;

		Random randInst = new Random();

		for (int n = 0; n < 100; n++) {

			double rand = randInst.nextGaussian();
			double val = mean + (sd * rand);
			if (val < seqlen) {
				val = 100;
			}
			System.out.println(val);

		}

	}

}
