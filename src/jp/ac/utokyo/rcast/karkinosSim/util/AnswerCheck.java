package jp.ac.utokyo.rcast.karkinosSim.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;

public class AnswerCheck {
	
	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();

		
		optionlist.add(getOption("answer", "answer", true,
				"answer vcf file", true));

		optionlist.add(getOption("check", "check", true, "vcf file to check", true));

		return optionlist;
	}
	


	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}


	public static void main(String[] arg) {
		
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
			help.printHelp("karkinosSim.jar checkAnswer", opts, true);
			return;
		}
		
		String answer = cl.getOptionValue("answer");
		String check = cl.getOptionValue("check");
		
		exec(answer, check);
		

//		 String[] ids = { "tp0", "tp5", "tp10", "tp20", "tp30", "tp40", "tp50",
//		 	"tp60", "tp70", "tp80", "tp90", "tp95", "tp100" };
//
//	//	 String[] ids ={"tp5"};
//		 String glf ="/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP3.vcf";
//		// Map<String, MutationBean> glfb = getMap(glf,-1);
//		 Map<String, MutationBean> glfb =null;
//		 
//		for (String id : ids) {
//
//			exec(id, glfb);
//		}

	}

	public static void exec(String answer,String check) {

		// String answer =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/THCC135_old/THCC_135T-THCC_135N_snvdata.vcf";
		// String check =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/result135/tp5/tp5/tp5_snvdata.vcf";
		//System.out.println(id);

		//String answer = "/GLUSTER_DIST/data/users/ueda/karkinosSim/sim1200/text1200_4.vcf";
		//String check = "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim1200_4recal/"
		//		+"/"+id+"/" + id + "/" + id + "_snvdata.vcf";

		//
//		 String check =
//		 "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim_1200_4mutectrecal/"+id+"/"
//		 +id+".vcf";

		// String check =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim1200_4varscantp/"+id+"_filter.snp";

		Map<String, MutationBean> ansM = getMap(answer, 0);
		Map<String, MutationBean> checkM = getMap(check, 0);

		Iterator<String> ite = ansM.keySet().iterator();
		int total = 0;
		int exsist = 0;
		int pass = 0;
		int pass2 = 0;
		int fnr = 0;
		while (ite.hasNext()) {

			String key = ite.next();
			MutationBean ans = ansM.get(key);
			if (ans.indel)
				continue;
//			if (glfb != null) {
//				if (glfb.containsKey(key)) {
//					continue;
//				}
//			}

			if (checkM.containsKey(key)) {
				exsist++;
				MutationBean mb = checkM.get(key);
				mb.setExsist(true);
				
				if (mb.pass) {
					pass++;
					mb.setPasscheck(true);
				}else {
					//System.out.println(mb.chr+"\t"+mb.pos+"\t"+mb.filter+"\t"+mb.filter2);
				}
				if (mb.pass2) {
					pass2++;
					mb.setPass2check(true);
				} 

			} else {
				// System.out.println("FN"+"\t"+ans.chr+"\t"+ans.pos+"\t"+ans.ref+"-"+ans.alt);
				 if(ans.ref.equals(ans.alt)){
					 fnr++;
				 }
			}
			total++;
		}
		System.out.println(total + "\t" + exsist + "\t" + pass + "\t" + pass2);
		//
		int passtotal = 0;
		int pass2total = 0;
		int passcheck = 0;
		int pass2check = 0;

		Iterator<String> ite2 = checkM.keySet().iterator();
		while (ite2.hasNext()) {

			String key = ite2.next();
			MutationBean checkb = checkM.get(key);
			if (checkb.indel)
				continue;
//			if (glfb != null) {
//				if (glfb.containsKey(key)) {
//					continue;
//				}
//			}

			if(!checkb.isExsist()){
				if(checkb.pass){
					//System.out.println(key + checkb.filter);
				}
			}
			
			if (checkb.pass) {
				passtotal++;
				if (checkb.passcheck) {
					passcheck++;
				}
			}else{
				if (checkb.isExsist()) {
					//System.out.println(key + checkb.filter);
				}
				
			}

			if (checkb.pass2) {
				pass2total++;
				if (checkb.pass2check) {
					pass2check++;
				}
			}

		}
		System.out.println(passtotal + "\t" + passcheck);
		System.out.println(pass2total + "\t" + pass2check);
		System.out.println("FNok="+fnr);
	}

	private static Map<String, MutationBean> getMap(String check, int mode) {
		Map<String, MutationBean> map = new LinkedHashMap<String, MutationBean>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(check));
			String str = br.readLine();
			while (str != null) {
				if (str.startsWith("#")) {
					str = br.readLine();
					continue;
				}
				if (str.startsWith("contig")) {
					str = br.readLine();
					continue;
				}
				if (str.startsWith("chrom")) {
					str = br.readLine();
					continue;
				}

				String[] stra = str.split("\t");
				String chr = stra[0].replaceAll("chr", "");
				String poskey = (chr + "-" + stra[1]);
				map.put(poskey, new MutationBean(stra, mode));
				str = br.readLine();
			}

			br.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return map;
	}

}
