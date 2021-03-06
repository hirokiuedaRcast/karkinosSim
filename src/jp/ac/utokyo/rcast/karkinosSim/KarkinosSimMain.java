package jp.ac.utokyo.rcast.karkinosSim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
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

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import jp.ac.utokyo.rcast.karkinosSim.refs.CapHolderSim;

public class KarkinosSimMain extends ReadWriteBase {

	public static void main(String[] arg) {

		try {
			exec(arg);
			//exec(arg);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static final int MODE_SNP = 1;
	public static final int MODE_Artifitial = 2;
	public static final int MODE_Artifitial_SUB = 3;

	public static void exec(String tp) throws Exception {

		String bed = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		String tb = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
//		 String nb1 =
//		 "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/THCC_135T-THCC_135N_normal_pre_0.bam";
//		 String nb2 =
//		 "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/THCC_135T-THCC_135N_normal_pre_1.bam";
		// String id ="test";
		// String outdir =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200_4/";
		 String SNPlist =
		 "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/tp100_normalsnp.vcf";
		// String artifitial =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/text1200_4.vcf";
		 String CNV
		 ="/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/testCNV.txt";

		String nb1 = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bamsim/testout150_1.bam";
		String nb2 = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bamsim/testout12_15_150_2.bam";
		String id = "test";
		String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim7/sim2000/";
		//String SNPlist = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/sim1200/THCC_12-15T-THCC_12-15N_normalsnp.vcf";
//		String artifitial = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/sim1200/text1200_n1.vcf";
//		String artifitial_sub = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/sim1200/text1200_n2_sp.vcf";
		
		String artifitial = "/GLUSTER_DIST/data/users/ueda/karkinosSim6/text1000.vcf";
		String artifitial_sub = "/GLUSTER_DIST/data/users/ueda/karkinosSim6/text1000_sub.vcf";
		
		//String CNV = "/GLUSTER_DIST/data/users/ueda/karkinosSim6/testCNV6.txt";

		//double tc[] = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1 };
		//double tcc=0.95;
		double tcc = Double.parseDouble(tp);
		// double tc[] = {0,0.95};
		//for (double tcc : tc) {

			id = "test" + tcc;
			String od = outdir + "/" + "tp_" + (int) (tcc * 100);
			File f = new File(od);
			if (!f.exists()) {
				f.mkdirs();
			}
			_exec(bed, tb, nb1, nb2, id, od, SNPlist, artifitial,artifitial_sub, CNV, tcc);

		//}
	}
	
	
	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();

		
		optionlist.add(getOption("id", "id", true,
				"id for this similation", true));
		
		optionlist.add(getOption("target", "target bed", true,
				"target bed file", true));
		
		optionlist.add(getOption("r", "ref", true,
				"2bit reference file", true));
		
		optionlist.add(getOption("SNP", "SNPVcf", true,
				"SNP list exclude SNP positions", true));
		
		optionlist.add(getOption("SNVgen", "SNVgen", true,
				"SNV to generate, vcf", true));
		
		optionlist.add(getOption("SNVgenSub", "SNVgenSub", true,
				"subpolulation SNV to generate, vcf", true));
		
			
		optionlist.add(getOption("CNVgen", "CNVgen", true,
				"CNV to generate, vcf", true));
		
		optionlist.add(getOption("tc", "tumorContent", true,
				"tumor content ration for similation", true));
		
		
				
		optionlist.add(getOption("normalbam1", "normalbam1", true,
				"normalbam1", true));
		

		optionlist.add(getOption("normalbam2", "normalbam2", true,
				"normalbam2", true));

		
		

		
		optionlist.add(getOption("out", "out", true, "output directory", true));
		
		

		return optionlist;
	}
	

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	
	public static void exec(String[] arg) throws Exception {
		
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
			help.printHelp("karkinosSim.jar similate", opts, true);
			return;
		}
		
		
		//optionlist.add(getOption("id", "id", true,
		String id = cl.getOptionValue("id");
		//optionlist.add(getOption("target", "target bed", true,
		String bed = cl.getOptionValue("target");					
		//optionlist.add(getOption("r", "ref", true,
		String tb = cl.getOptionValue("r");	
		//optionlist.add(getOption("SNP", "SNPVcf", true,
		String SNPlist = cl.getOptionValue("SNP");	
		//optionlist.add(getOption("SNVgen", "SNVgen", true,
		String artifitial = cl.getOptionValue("SNVgen");	
		//optionlist.add(getOption("SNVgenSub", "SNVgenSub", true,
		String artifitial_sub = cl.getOptionValue("SNVgenSub");
		//optionlist.add(getOption("CNVgen", "CNVgen", true,
		String CNV 	= cl.getOptionValue("CNVgen");	
		//optionlist.add(getOption("tc", "tumorContent", true,
					
		
				
		
		//optionlist.add(getOption("normalbam1", "normalbam1", true,
		String nb1 = 	cl.getOptionValue("normalbam1");	

		//optionlist.add(getOption("normalbam2", "normalbam2", true,
		String nb2 = 	cl.getOptionValue("normalbam2");
		
						
		//optionlist.add(getOption("out", "out", true, "output directory", true));
		
		String outdir = 	cl.getOptionValue("out");
	

		double tcc = Double.parseDouble(cl.getOptionValue("tc"));
		
		_exec(bed, tb, nb1, nb2, id, outdir, SNPlist, artifitial,artifitial_sub, CNV, tcc);
		
	}

	private static void _exec(String bed, String tb, String nb1, String nb2,
			String id, String outdir, String SNPlist, String artifitial,String artifitial_sub,
			String CNV, double tc) throws IOException, Exception {
		//
		int librarysizemax = 1000;
		//
		CapHolderSim capholder = new CapHolderSim();
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(tb));
		CNVHolder cnv = new CNVHolder(new File(CNV));
		// reads target region
		capholder.loadTargetBed(bed, tgr);
		Map<String, SNPHolder> snpHolder = new LinkedHashMap<String, SNPHolder>();

		// vcf
		regSNP(snpHolder, capholder, SNPlist, MODE_SNP);
		// check AF;

		// reg Artifitial mutation
		regSNP(snpHolder, capholder, artifitial, MODE_Artifitial);
		
		regSNP(snpHolder, capholder, artifitial_sub, MODE_Artifitial_SUB);

		Iterator<String> ite = snpHolder.keySet().iterator();
		while(ite.hasNext()){
			
			String k = ite.next();
			System.out.println(k);
			if(k.contains("10")){
				
				SNPHolder sh = snpHolder.get(k);
				for(int n:sh.snpmap.keySet()){
					
					System.out.println(n);
					if(n==93726528){
						Object o = sh.snpmap.get(n);
						
						System.out.println("here");
						
					}
				}
				
			}
			
		}
		
		// check possible chain
		int sumlink = checkChain(snpHolder, librarysizemax);

		// analyze alle get allel A,Allele B map
		analyseAllele(snpHolder, nb1, nb2);

		//
		// Similator.similate(snpHolder, nb1, nb2,tgr,
		// capholder,cnv,tc,librarysizemax,id,outdir);

		Similator.similate(snpHolder, nb1, nb2, tgr, capholder, cnv, tc,
				librarysizemax, id, outdir);

		System.out.println(sumlink);
	}

	private static void analyseAllele(Map<String, SNPHolder> map, String nb1,
			String nb2) {
		//
		Iterator<String> ite = map.keySet().iterator();
		//
		SAMFileReader normalbam1 = getReader(nb1);
		SAMFileReader normalbam2 = getReader(nb2);

		while (ite.hasNext()) {
			String chr = ite.next();
			System.out.println(chr);
			SNPHolder snpHolder = map.get(chr);
			Map<Integer, AlleleFragment> afmap = snpHolder.getAfmap();
			Iterator<Integer> iteaf = afmap.keySet().iterator();
			while (iteaf.hasNext()) {

				int key = iteaf.next();
				AlleleFragment afg = afmap.get(key);
				Map<Integer, AllelleBase> allelemap = new LinkedHashMap<Integer, AllelleBase>();
				//
				LinkedHashMap<String, SamBean> readsmap = new LinkedHashMap<String, SamBean>();
				addSamList(chr, afg, normalbam1, readsmap);
				addSamList(chr, afg, normalbam2, readsmap);
				//
				alleleSeek(chr, afg, readsmap, allelemap);
				afg.setAllelemap(allelemap);
				// System.out.println(afg.getALMapStr());

			}

		}

		normalbam1.close();
		normalbam2.close();

	}

	private static void alleleSeek(String chr, AlleleFragment afg,
			LinkedHashMap<String, SamBean> readsmap,
			Map<Integer, AllelleBase> allelemap) {

		SNPInfoBean prev = null;

		if (afg.snpList.size() == 1) {

			//
			SNPInfoBean info = afg.snpList.get(0);
			AllelleBase ab1 = new AllelleBase();
			ab1.alleleA = info.ref;
			ab1.alleleB = info.alt;
			allelemap.put(info.pos, ab1);
			return;

		}

		for (SNPInfoBean info : afg.snpList) {

			if (info.modeSnp == 2) {
				continue;
			}
			//
			if (prev == null) {
				// no nothing

			} else {

				assign(chr, prev, info, readsmap, allelemap);

			}
			prev = info;
		}

	}

	private static void assign(String chr, SNPInfoBean prev, SNPInfoBean info,
			LinkedHashMap<String, SamBean> readsmap,
			Map<Integer, AllelleBase> allelemap) {

		//
		//
		int pos1 = prev.pos;
		boolean indel1 = prev.isIndel();
		int pos2 = info.pos;
		boolean indel2 = info.isIndel();
		//
		List<SamBean> listContain = contain(pos1, pos2, readsmap);
		//
		List<Counter> list = new ArrayList<Counter>();

		for (SamBean sb : listContain) {

			IndelInfo ii1 = new IndelInfo();
			char nuc1 = getNuc(pos1, sb.sam1, ii1);
			if (nuc1 == 'N') {
				nuc1 = getNuc(pos1, sb.sam2, ii1);
			}
			String s1 = nuc1 + "";
			if (indel1) {
				//
				s1 = ii1.indel ? "Indel" : "ref";

			}

			IndelInfo ii2 = new IndelInfo();
			char nuc2 = getNuc(pos2, sb.sam1, ii2);
			if (nuc2 == 'N') {
				nuc2 = getNuc(pos2, sb.sam2, ii2);
			}
			String s2 = nuc2 + "";
			if (indel2) {
				//
				s2 = ii2.indel ? "Indel" : "ref";
			}
			String key = s1 + "/" + s2;
			//
			Counter cnt = search(list, key);
			if (cnt == null) {
				cnt = new Counter(key);
				list.add(cnt);
			} else {
				cnt.inc();
			}

		}
		List<Counter> topTwo = getTopTwo(list);
		// debug
		// for(Counter ct:topTwo){
		// System.out.println(ct.key+"|"+ct.n);
		// }

		reg(prev, info, topTwo, allelemap);

	}

	private static void reg(SNPInfoBean prev, SNPInfoBean info,
			List<Counter> topTwo, Map<Integer, AllelleBase> allelemap) {

		//
		int pos1 = prev.pos;
		int pos2 = info.pos;
		Map<String, String> map = toMap(topTwo);

		//
		if (!allelemap.containsKey(pos1)) {
			// init
			AllelleBase ab1 = new AllelleBase();
			AllelleBase ab2 = new AllelleBase();
			Iterator<String> ite = map.keySet().iterator();
			int n = 0;
			while (ite.hasNext()) {
				String s1 = ite.next();
				String val = map.get(s1);
				if (n == 0) {
					ab1.alleleA = s1;
					ab2.alleleA = val;
				} else {
					ab1.alleleB = s1;
					ab2.alleleB = val;
				}
				n++;
			}
			allelemap.put(pos1, ab1);
			if (ab1.alleleA == null) {
				ab1.alleleA = prev.ref;
			}
			if (ab1.alleleB == null) {
				ab1.alleleB = prev.alt;
			}

			if (ab2.alleleA == null) {
				ab2.alleleA = info.ref;
			}
			if (ab2.alleleB == null) {
				ab2.alleleB = info.alt;
			}
			allelemap.put(pos2, ab2);
		} else {

			AllelleBase ab1 = allelemap.get(pos1);
			AllelleBase ab2 = new AllelleBase();

			ab2.alleleA = map.get(ab1.alleleA);
			ab2.alleleB = map.get(ab1.alleleB);
			if (ab2.alleleA == null) {
				ab2.alleleA = info.ref;
			}
			if (ab2.alleleB == null) {
				ab2.alleleB = info.alt;
			}
			allelemap.put(pos2, ab2);

		}

	}

	private static Map<String, String> toMap(List<Counter> topTwo) {

		Map<String, String> map = new HashMap<String, String>();
		for (Counter ct : topTwo) {

			String key1 = ct.key;
			String s1 = key1.substring(0, key1.indexOf("/"));
			String s2 = key1.substring(key1.indexOf("/") + 1);
			map.put(s1, s2);

		}
		return map;
	}

	private static List<Counter> getTopTwo(List<Counter> list) {
		List<Counter> l = new ArrayList<Counter>();
		Collections.sort(list, new CounterIterator());
		for (int n = 0; n < 2 && n < list.size(); n++) {
			l.add(list.get(n));
		}
		return l;
	}

	private static Counter search(List<Counter> list, String key) {

		for (Counter cnt : list) {
			if (cnt.key.equals(key)) {
				return cnt;
			}
		}
		return null;
	}

	private static char getNuc(int pos, SAMRecord sam, IndelInfo ii) {

		if (sam == null)
			return 'N';
		int mutationidx = KarkinosSimUtils.getCharIdx(pos, sam, ii);
		int readslen = sam.getReadLength();
		if ((mutationidx >= 0) && (mutationidx < readslen)) {
			char ch = (char) sam.getReadBases()[mutationidx];
			return ch;
			//
		}
		return 'N';
	}

	private static List<SamBean> contain(int pos1, int pos2,
			LinkedHashMap<String, SamBean> map) {

		//
		Set<Entry<String, SamBean>> set = map.entrySet();
		List<SamBean> list = new ArrayList<SamBean>();

		for (Entry<String, SamBean> et : set) {
			SamBean sb = et.getValue();
			boolean contain = contain(sb, pos1, pos2);
			if (contain) {
				list.add(sb);
			}
		}
		return list;
	}

	private static boolean contain(SamBean sb, int pos1, int pos2) {

		//
		boolean pos1contain = sb.contain(pos1);
		boolean pos2contain = sb.contain(pos2);
		return pos1contain && pos2contain;
	}

	private static void addSamList(String chr, AlleleFragment afg,
			SAMFileReader normalbam, LinkedHashMap<String, SamBean> map) {

		CloseableIterator<SAMRecord> normarIte = normalbam.query(chr,
				afg.start, afg.end, false);

		while (normarIte.hasNext()) {

			SAMRecord sam = normarIte.next();
			String readname = sam.getReadName();

			SamBean sb = null;
			if (map.containsKey(readname)) {
				sb = map.get(readname);
			} else {
				sb = new SamBean();
				map.put(readname, sb);
			}
			//
			//
			sb.add(sam);

		}
		normarIte.close();

	}

	private static int checkChain(Map<String, SNPHolder> map, int librarysizemax) {

		Iterator<String> ite = map.keySet().iterator();
		int totalfrag = 0;
		while (ite.hasNext()) {

			//
			String chr = ite.next();
			// System.out.println("chr="+chr);
			SNPHolder snpHolder = map.get(chr);
			snpHolder.checkLink(librarysizemax);

			// debug
			// Map<Integer,AlleleFragment> afmap = snpHolder.getAfmap();
			// Iterator<Integer> iteaf = afmap.keySet().iterator();
			//
			// while(iteaf.hasNext()){
			//
			// int key = iteaf.next();
			// AlleleFragment afg = afmap.get(key);
			// System.out.print(chr+"\t"+afg.start+"\t"+afg.end+"\t"+(afg.end-afg.start)+"\t"+afg.count);
			// System.out.println("\t"+afg.getAInfoStr());
			// totalfrag++;
			// }

		}
		return totalfrag;

	}

	private static void regSNP(Map<String, SNPHolder> snpHolder,
			CapHolderSim capholder, String dbSNP, int modeSnp)
			throws IOException {

		int margin = 100;
		FileInputStream fis = new FileInputStream(dbSNP);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		int idxchr = 0;
		int idxpos = 1;

		int cnt = 0;
		try {

			int totalcnt = 0;
			int init = fis.available();
			String chr = "";
			for (;;) {
				int point = init - fis.available();
				String line = br.readLine();

				if (line == null)
					break;
				if (line.startsWith("#")) {
					continue;
				}
				totalcnt++;
				String[] sa = line.split("\t");
				String _chr = sa[idxchr];
				if (!_chr.contains("chr")) {
					_chr = "chr" + _chr;
				}
				if (!chr.equals(_chr)) {
					System.out.println(_chr + "\t" + totalcnt + "\t" + cnt);
				}
				chr = _chr;
				int pos = Integer.parseInt(sa[idxpos]);

				float af = 0.5f;
				String id = sa[2];
				String ref = sa[3];
				String alt = sa[4];

				if (sa.length > 7) {
					String is = sa[7];
					if (is.contains("AF")) {
						String afs = is.substring(is.indexOf("AF=") + 3,
								is.indexOf("AF=") + 6);
						// System.out.println(afs);
						try {
							af = Float.parseFloat(afs);
						} catch (Exception ex) {

						}
					}
				}

				boolean capr = (capholder.getOverlapping(chr, pos - margin, pos
						+ margin) != null);
				boolean reg = capr && af <= 0.8;
				if (modeSnp == 2) {
					reg = true;
				}
				if (reg) {
					// reg
					SNPHolder sh = snpHolder.get(chr);
					if (sh == null) {
						sh = new SNPHolder();
						snpHolder.put(chr, sh);
					}//

					sh.add(chr, pos, af, id, ref, alt, modeSnp);
					cnt++;
				}
				// if(cnt>100)break;
			}

		} finally {
			br.close();
		}
		System.out.println(cnt);

	}

}
