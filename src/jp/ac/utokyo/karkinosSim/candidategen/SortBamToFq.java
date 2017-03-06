package jp.ac.utokyo.karkinosSim.candidategen;

import java.util.LinkedHashMap;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinosSim.FastQwriter;
import jp.ac.utokyo.rcast.karkinosSim.SamBean;

public class SortBamToFq extends ReadWriteBase {
	
//	final static String inbam = "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/THCC_135T-THCC_135N_normal_pre_2.bam";
//	final static String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200_2/ref/";
	
	final static String inbam = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bamsim/testout12_15_0.bam";
	final static String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/sim1200/ref";
	
	final static String id="ref"; 
	
	public static void main(String[] arg){
		
		SortBamToFq inst = new SortBamToFq();
		try {
			inst.exec();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void exec() throws Exception {
		
		//
		SAMFileReader normalbam = getReader(inbam);
		FastQwriter fqw = new FastQwriter(outdir + "/fq/", id, 4000000);
		//
		CloseableIterator<SAMRecord> ite = normalbam.iterator();
		LinkedHashMap<String, SamBean> readsmap = new LinkedHashMap<String, SamBean>();
		while(ite.hasNext()){
			
			
			SAMRecord sr = ite.next();				
			SamBean sb = add(readsmap, sr);
			if(sb.bothExsist()){
				fqw.add(sb);
				readsmap.remove(sr.getReadName());
			}
			
		}		
		ite.close();
		normalbam.close();
		fqw.close();
	}
	
	
	private static SamBean add(LinkedHashMap<String, SamBean> map, SAMRecord sr) {

		String readname = sr.getReadName();
		SamBean sb = null;
		if (map.containsKey(readname)) {
			sb = map.get(readname);
		} else {
			sb = new SamBean();
			map.put(readname, sb);
		}
		//
		//
		sb.add(sr);
		return sb;

	}
}
