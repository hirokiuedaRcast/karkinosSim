package jp.ac.utokyo.karkinosSim.candidategen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import jp.ac.utokyo.rcast.karkinosSim.CNVHolder;
import jp.ac.utokyo.rcast.karkinosSim.refs.CapHolderSim;

public class RomdomAssign {
	
	
	public static void main(String[] arg){
		
		//
		try {
			
			exec(arg);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();

		
		optionlist.add(getOption("target", "target bed", true,
				"target bed file", true));
		
		optionlist.add(getOption("r", "ref", true,
				"2bit reference file", true));
		
		optionlist.add(getOption("lowdepth", "lowdepth", true,
				"low depth vcf", true));

		
		optionlist.add(getOption("outvcf", "outvcf", true, "output vcf file", true));

		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	public static void exec(String[] arg) throws IOException {

		//
//		String bed = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
//		String tb = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
//		String outvcf = "/GLUSTER_DIST/data/users/ueda/karkinosSim5/text1000_sub.vcf";
//		String lowdap1 ="/GLUSTER_DIST/data/users/ueda/karkinosSim/result135/tp10/tp10/tp10_tumor_lowcov.bed";
//		//String lowdap2="/GLUSTER_DIST/data/users/ueda/karkinosSim/result135/tp10/tp10/tp10_normal_lowcov.bed";
		
		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOption();
		Options opts = new Options();
		for (Option opt : optionList) {
			opts.addOption(opt);
		}

		CommandLine cl = null;
		try {
			cl = parcer.parse(opts, arg);
		} catch (ParseException e1) {
			System.out.println(e1.getMessage());
			HelpFormatter help = new HelpFormatter();
			help.setOptionComparator(new OptionComparator(optionList));
			help.printHelp("karkinos.jar varidate", opts, true);
			return;
		}
		
		
		String bed = cl.getOptionValue("target");
		String tb = cl.getOptionValue("r");
		String outvcf = cl.getOptionValue("outvcf");
		String lowdap1 = cl.getOptionValue("lowdepth");
		
		Map<String,Set<Integer>> map = new HashMap<String,Set<Integer>>();
		
		load(map,lowdap1);
		//load(map,lowdap2);
		
		
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
				new FileOutputStream(outvcf)));
		writeHeader(bw);
		
		int numSNV = 1000;
		int numins = 0;
		int numdel = 0;
		
		float indelfrac = 0.15f;

		int maxIndelLength = 20;

		int mutationmindist = 100;
		int librarysizemax = 1000;

		//
		CapHolderSim capholder = new CapHolderSim();
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(tb));
		// reads target region
		try {
			capholder.loadTargetBed(bed, tgr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// /
		TreeMap<Long, CapInterval> indexmap = new TreeMap<Long, CapInterval>();
		long total = setIndexMap(indexmap, capholder);
		
		List<Position> list = new ArrayList<Position>();
		
		for (int n = 0; n < numSNV; n++) {

			Position pos = random(indexmap, total);
			if(!depthcheck(pos,map)){
				list.add(pos);
			}else{
				n--;
			}
		}	
		for (int n = 0; n < numins; n++) {

			Position pos = random(indexmap, total);
			pos.setIndel(true);
			pos.setInsersion(true);
			if(!depthcheck(pos,map)){
			  list.add(pos);
			}else{
			  n--;	
			}
		}	
		for (int n = 0; n < numins; n++) {

			Position pos = random(indexmap, total);
			pos.setIndel(true);
			if(!depthcheck(pos,map)){
			 list.add(pos);
			}else{
				n--;
			}
		}	
		
		sort(list);
		int num = 0;
		for (Position pos:list) {

			num++;
			String ref = "";
			String alt = "";
			//
			
			
			// //
			int len = 0;
			if (pos.isIndel()) {
				//
				len = getLength(maxIndelLength);
			}
			ref = tgr.getGenomicSeq(pos.chr, pos.pos, pos.pos + (len), true);

			if (pos.isIndel()&&pos.isInsersion()) {
				//
				ref = tgr.getGenomicSeq(pos.chr, pos.pos, pos.pos, true);
				alt = tgr.getGenomicSeq(pos.chr, pos.pos, pos.pos, true);
				alt = alt + randstr(len);
				
			} else if (pos.isIndel()) {
				//del
				
				alt = tgr.getGenomicSeq(pos.chr, pos.pos, pos.pos, true);

			} else {
				//
				ref = tgr.getGenomeNuc(pos.chr, pos.pos, true) + "";
				alt = randDiffstr(ref);
				
			}
			ref = ref.toUpperCase();
			alt = alt.toUpperCase();
			
			System.out.println(pos.chr+"\t"+ pos.pos+"\t"+ref+"\t"+alt);
			bw.append(pos.chr.replace("chr","")+"\t"+
					pos.pos+"\t"+
					"sum"+num+"\t"+
					ref+"\t"+
					alt+"\t"+
					"100"+"\t"+
					"PASS"+"\t"+
					"noinfo"+"\t"+
					"PASS"+"\t"+
					"\n");
			
		}
		bw.close();

	}

	private static void load(Map<String, Set<Integer>> map, String lowdap1) {
		
		//
		if(map==null){
			map = new HashMap<String,Set<Integer>>();
		}
		//
		FileInputStream fis;
		try {
			
			fis = new FileInputStream(lowdap1);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
			String line = null;
			while( (line = br.readLine())!=null ){
				
				if(line.startsWith("#"))continue;
				String[] sa = line.split("\t");
				String chr = sa[0];
				int start = Integer.parseInt(sa[1]);
				int end = Integer.parseInt(sa[2]);
				Set<Integer> s= null;
				if(map.containsKey(chr)){
					s = map.get(chr);
				}else{
					s = new HashSet<Integer>();
					map.put(chr, s);
				}
				//
				for(int n= start;n<=end;n++){
					s.add(n);
				}
				
			}
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
		
	}

	private static boolean depthcheck(Position pos,
			Map<String, Set<Integer>> map) {
		
		String chr = pos.chr;
		int posn = pos.pos;
		
		if(map.containsKey(chr)){
			
			Set s = map.get(chr);
			if(s.contains(posn)){
				return true;
			}
			
		}
		return false;
	
	}

	private static void writeHeader(BufferedWriter bw) throws IOException {
		// TODO Auto-generated method stub
		bw.append("##fileformat=VCFv4.1"+ "\n");
		bw.append("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FILTER2	Normal	Tumor	"+ "\n");
		
	}

	private static void sort(List<Position> list) {
		
		Collections.sort(list, new MyComparator());
		
	}

	static class MyComparator implements Comparator<Position>{

		public int compare(Position p0, Position p1) {
			
			String chr0 = p0.getChr().replace("chr","");
			String chr1 = p1.getChr().replace("chr","");
			
			if(chr0.equals(chr1)){
				return p0.pos < p1.pos ? -1:1;
			}else{
				return chr0.compareTo(chr1);
			}			
		}
		
		
	}
	
	private static String randDiffstr(String ref) {

		
		String nuc = randNuc()+"";
		while(nuc.equals(ref)){
			nuc = randNuc()+"";			
		}
		return nuc;

	}

	private static String randstr(int len) {

		StringBuffer sb = new StringBuffer();
		for (int n = 1; n <= len; n++) {

			sb.append(randNuc());

		}
		return sb.toString();

	}

	private static char randNuc() {

		double rand = Math.random();
		double unit = 0.25;
		if (rand <= unit) {
			return 'A';
		} else if (rand <= (unit * 2)) {
			return 'T';
		} else if (rand <= (unit * 3)) {
			return 'C';
		} else {
			return 'G';
		}
	}

	private static int getLength(int maxIndelLength) {

		//
		double rand = Math.random() * maxIndelLength;
		int len = (int) rand;
		if (len == 0) {
			len = 1;
		}
		return len;
	}

	private static Position random(TreeMap<Long, CapInterval> indexmap, long total) {

		long randompos = Math.round(total * Math.random());
		long floorkey = indexmap.floorKey(randompos);
		CapInterval cp = indexmap.get(floorkey);
		Position pos = new Position();
		pos.chr = cp.getChr();
		pos.pos = (int) (cp.getStart() + (randompos - floorkey));
		return pos;
	}

	private static long setIndexMap(TreeMap<Long, CapInterval> indexmap,
			CapHolderSim capholder) {

		long total = 1;
		Map<String, TreeMap<Integer, CapInterval>> map = capholder.getMap();
		Iterator<String> ite = map.keySet().iterator();
		while (ite.hasNext()) {

			//
			String chr = ite.next();
			TreeMap<Integer, CapInterval> cim = map.get(chr);
			Set<Entry<Integer, CapInterval>> set = cim.entrySet();
			for (Entry<Integer, CapInterval> et : set) {
				//
				int intervallen = et.getValue().getLength();
				indexmap.put(total, et.getValue());
				total = total + intervallen;
			}

		}
		return total;
	}

}