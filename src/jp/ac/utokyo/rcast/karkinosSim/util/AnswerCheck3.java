package jp.ac.utokyo.rcast.karkinosSim.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinosSim.Counter;

public class AnswerCheck3 {

	public static void main(String[] sa) {

		 String[] ids = { "tp0", "tp05", "tp10", "tp20","tp30",  "tp40", "tp50",
		 "tp60",
		 "tp70", "tp80", "tp90", "tp95", "tp100" };

//		String[] ids = { "tp05" };

		String glf = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP4.vcf";
		Map<String, MutationBean> glfb = getMap(glf, -1);

		// Map<String, MutationBean> glfb =null;

		String answer = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/sim1200/text1200_n1.vcf";
		Map<String, MutationBean> ansM = getMap(answer, 0);

		String artifitial_sub = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/sim1200/text1200_n2_sp.vcf";

		Map<String, MutationBean> ansM2 = getMap(artifitial_sub, 0);

		Iterator<String> ite = ansM.keySet().iterator();
		List<String> removekeys = new ArrayList<String>();

		while (ite.hasNext()) {
			String key = ite.next();
			boolean remove = false;
			if (glfb.containsKey(key)) {
				remove = true;
			}
			MutationBean ans = ansM.get(key);
			if (ans.ref == ans.alt) {
				remove = true;
			}
			if (remove) {
				removekeys.add(key);
			}

		}
		for (String key : removekeys) {
			ansM.remove(key);
		}

		Iterator<String> ite2 = ansM2.keySet().iterator();
		List<String> removekeys2 = new ArrayList<String>();

		while (ite2.hasNext()) {
			String key = ite2.next();
			boolean remove = false;
			if (glfb.containsKey(key)) {
				remove = true;
			}
			MutationBean ans = ansM2.get(key);
			if (ans.ref == ans.alt) {
				remove = true;
			}
			if (remove) {
				removekeys2.add(key);
			}

		}
		for (String key : removekeys2) {
			ansM2.remove(key);
		}

		// remove normal panel

		Set<String> s = new HashSet<String>();
		Map<String, AnswerBean> ab = new LinkedHashMap<String, AnswerBean>();
		System.out.println(ansM.size());

		for (String id : ids) {

			String check = "/GLUSTER_DIST/data/users/ueda/karkinosSim4/Mutect/"
					+ id + "/" + id + ".vcf";

			// String check =
			// "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim1200_4varscantp/"+id+"_filter.snp";

			Map<String, MutationBean> checkM = getMap(check, 1);
			Iterator<String> itec = checkM.keySet().iterator();
			List<String> remlist = new ArrayList<String>();
			while (itec.hasNext()) {
				String key = itec.next();
				if (glfb.containsKey(key)) {
					remlist.add(key);
				}
			}
			for (String ss : remlist) {
				checkM.remove(ss);
			}

			exec(id, checkM, glfb, ansM, ansM2, s, ab);

			String check2 = "/GLUSTER_DIST/data/users/ueda/karkinosSim4/karkinos/"
					+ id + "/" + id + "/" + id + "_snvdata.vcf";

			// String check =
			// "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim1200_4varscantp/"+id+"_filter.snp";

			Map<String, MutationBean> checkM2 = getMap(check2, 0);
			Iterator<String> itec2 = checkM2.keySet().iterator();
			List<String> remlist2 = new ArrayList<String>();
			while (itec2.hasNext()) {
				String key = itec2.next();
				if (glfb.containsKey(key)) {
					remlist2.add(key);
				}
			}
			for (String ss : remlist2) {
				checkM2.remove(ss);
			}

			exec(id, checkM2, glfb, ansM, ansM2, s, ab);

		}
		// display
		Iterator<String> ite3 = ab.keySet().iterator();
		while (ite3.hasNext()) {

			String id = ite3.next();
			AnswerBean abl = ab.get(id);
			//
			for (String ss : s) {

				Counter c1 = abl.amc.get(ss);
				Counter c2 = abl.mc.get(ss);
				int n = 0;
				int m = 0;
				if (c1 != null) {
					m = c1.getN();
				}
				if (c2 != null) {
					n = c2.getN();
				}

				float fp = (float) ((double) n / (double) m);
				System.out.println(id + "\t" + ss + "\t" + m + "\t" + fp);

			}

		}

	}

	public static void exec(String id, Map<String, MutationBean> checkM,
			Map<String, MutationBean> glfb, Map<String, MutationBean> ansM,
			Map<String, MutationBean> ansM2, Set<String> set,
			Map<String, AnswerBean> ab) {

		// String answer =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/THCC135_old/THCC_135T-THCC_135N_snvdata.vcf";
		// String check =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/result135/tp5/tp5/tp5_snvdata.vcf";
		System.out.println(id);

		// String check =
		// "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim1200_4recal/"
		// + "/" + id + "/" + id + "/" + id + "_snvdata.vcf";

		//
		Set<String> hit = new HashSet<String>();
		Set<String> ansindel = new HashSet<String>();

		int exsist = 0;
		int exsistf = 0;
		int miss = 0;
		Iterator<String> ite = ansM.keySet().iterator();
		while (ite.hasNext()) {

			String key = ite.next();
			if (ansM.get(key).isIndel()) {
				ansindel.add(key);
				continue;
			}

			if (checkM.containsKey(key)) {

				if (checkM.get(key).pass2) {
					exsist++;
				} else {
					exsistf++;
//					System.out.println("Main p " + checkM.get(key).chr + "\t"
//							+ checkM.get(key).pos + "\t"
//							+ checkM.get(key).filter2);
				}
				
//				if ((checkM.get(key).pass) && (!checkM.get(key).pass2)){
//					System.out.println("Main p " + checkM.get(key).chr + "\t"
//					+ checkM.get(key).pos + "\t"
//					+ checkM.get(key).filter2);
//				}
				
				hit.add(key);
			} else {
				miss++;
			}

		}
		int exsist2 = 0;
		int exsistf2 = 0;
		int miss2 = 0;
		Iterator<String> ite2 = ansM2.keySet().iterator();
		while (ite2.hasNext()) {

			String key = ite2.next();
			if (ansM2.get(key).isIndel()) {
				ansindel.add(key);
				continue;
			}
			if (checkM.containsKey(key)) {
				if (checkM.get(key).pass2) {
					exsist2++;
				} else {
					exsistf2++;
//					System.out.println("Sub p " + checkM.get(key).chr + "\t"
//							+ checkM.get(key).pos + "\t"
//							+ checkM.get(key).filter2);
				}
				
//				if ((checkM.get(key).pass) && (!checkM.get(key).pass2)){
//					System.out.println("Main p " + checkM.get(key).chr + "\t"
//					+ checkM.get(key).pos + "\t"
//					+ checkM.get(key).filter2+ "\t");
//				}
				hit.add(key);
			} else {
				miss2++;
			}

		}
		System.out.print("e1=" + exsist + "\t" + "ef=" + exsistf + "\t"
				+ "miss=" + miss);
		System.out.println("\t e2=" + exsist2 + "\t" + "e2f=" + exsistf2 + "\t"
				+ "miss=" + miss2);

		int total = 0;
		Iterator<String> ite3 = checkM.keySet().iterator();
		while (ite3.hasNext()) {

			String key = ite3.next();
			if (checkM.get(key).isIndel()) {
				continue;
			}
			if (checkM.get(key).pass2) {
				if(!ansindel.contains(key)){
					total++;
					if (!hit.contains(key)) {
						
//						if(key.contains("397352")){
//							System.out.println(key);
//						}
						System.out.println("fp " + checkM.get(key).chr + "\t"
								+ checkM.get(key).pos + "\t"
								+ checkM.get(key).filter2+"\t"+ansindel.contains(key));
					}	
				}
				
			}			
			
		}
		System.out.println("total=" + total);
		// ans
		// Map<String,Counter> mc = new TreeMap<String,Counter>();
		// Map<String,Counter> mcu = new TreeMap<String,Counter>();
		//
		// //All
		// Map<String,Counter> amc = new TreeMap<String,Counter>();
		// Map<String,Counter> amcu = new TreeMap<String,Counter>();
		//
		// Iterator<String> ite = ansM.keySet().iterator();
		// int total = 0;
		// int exsist = 0;
		// int pass = 0;
		// int pass2 = 0;
		// int fnr = 0;
		// while (ite.hasNext()) {
		//
		// String key = ite.next();
		// MutationBean ans = ansM.get(key);
		// if (ans.indel)
		// continue;
		// if (glfb != null) {
		// if (glfb.containsKey(key)) {
		// continue;
		// }
		// }
		//
		// if (checkM.containsKey(key)) {
		// exsist++;
		// MutationBean mb = checkM.get(key);
		// mb.setExsist(true);
		//
		// if (mb.pass) {
		// pass++;
		// mb.setPasscheck(true);
		// } else {
		// //
		// System.out.println(mb.chr+"\t"+mb.pos+"\t"+mb.filter+"\t"+mb.filter2);
		// String s = mb.filter;
		// if(s.contains(",")){
		//
		// String[] sa = s.split(",");
		// for(String ss:sa){
		//
		// //
		// set.add(ss);
		// if(ss.length()>1){
		// Counter cnt = null;
		// if(mc.containsKey(ss)){
		// cnt = mc.get(ss);
		// }else{
		// cnt = new Counter();
		// mc.put(ss, cnt);
		// }
		// cnt.inc();
		// }
		//
		// }
		//
		// }else{
		// Counter cnt = null;
		// if(mcu.containsKey(s)){
		// cnt = mcu.get(s);
		// }else{
		// cnt = new Counter();
		// mcu.put(s, cnt);
		// }
		// cnt.inc();
		//
		// Counter cnt2 = null;
		// if(mc.containsKey(s)){
		// cnt2 = mc.get(s);
		// }else{
		// cnt2 = new Counter();
		// mc.put(s, cnt2);
		// }
		// cnt2.inc();
		// }
		//
		// }
		// if (mb.pass2) {
		// pass2++;
		// mb.setPass2check(true);
		// }
		//
		// } else {
		// //
		// System.out.println("FN"+"\t"+ans.chr+"\t"+ans.pos+"\t"+ans.ref+"-"+ans.alt);
		// if (ans.ref.equals(ans.alt)) {
		// fnr++;
		// }
		// }
		// total++;
		// }
		// System.out.println(total + "\t" + exsist + "\t" + pass + "\t" +
		// pass2);
		//
		//
		// //
		// int passtotal = 0;
		// int pass2total = 0;
		// int passcheck = 0;
		// int pass2check = 0;
		//
		// Iterator<String> ite2 = checkM.keySet().iterator();
		// int counttotal = 0;
		// int filterout = 0;
		// while (ite2.hasNext()) {
		//
		// String key = ite2.next();
		// MutationBean checkb = checkM.get(key);
		// if (checkb.indel)
		// continue;
		// if (glfb != null) {
		// if (glfb.containsKey(key)) {
		// continue;
		// }
		// }
		// counttotal++;
		//
		// if (!checkb.isExsist()) {
		// if (checkb.pass) {
		// // System.out.println(key);
		// }
		// }
		//
		// if (checkb.pass) {
		// passtotal++;
		// if (checkb.passcheck) {
		// passcheck++;
		// }
		// }else{
		// filterout++;
		// String s = checkb.filter;
		// if(s.contains("illumina")){
		// //System.out.println(s);
		// }
		// if(s.contains(",")){
		//
		// String[] sa = s.split(",");
		// for(String ss:sa){
		//
		// //
		// set.add(ss);
		// if(ss.length()>1){
		// Counter cnt0 = null;
		// if(amc.containsKey(ss)){
		// cnt0 = amc.get(ss);
		// }else{
		// cnt0 = new Counter();
		// amc.put(ss, cnt0);
		// }
		// cnt0.inc();
		// }
		//
		// }
		//
		// }else{
		// Counter cnt = null;
		// if(amcu.containsKey(s)){
		// cnt = amcu.get(s);
		// }else{
		// cnt = new Counter();
		// amcu.put(s, cnt);
		// }
		// cnt.inc();
		//
		// Counter cnt2 = null;
		// if(amc.containsKey(s)){
		// cnt2 = amc.get(s);
		// }else{
		// cnt2 = new Counter();
		// amc.put(s, cnt2);
		// }
		// cnt2.inc();
		// }
		// }
		//
		// if (checkb.pass2) {
		// pass2total++;
		// if (checkb.pass2check) {
		// pass2check++;
		// }
		// }
		//
		// }
		// System.out.println("initcand="+ counttotal+ "\t filter=" +
		// filterout);
		// System.out.println(passtotal + "\t" + passcheck);
		// System.out.println(pass2total + "\t" + pass2check);
		// System.out.println("FNok=" + fnr);
		// AnswerBean abean = new AnswerBean();
		// abean.setAmc(amc);
		// abean.setMc(mc);
		// ab.put(id, abean);
		//
		// // Iterator<String> itekey = amc.keySet().iterator();
		// // while(itekey.hasNext()){
		// //
		// // String k = itekey.next();
		// // int au = 0;
		// // try{
		// //
		// // au =amcu.get(k).getN();
		// //
		// // }catch(NullPointerException ex){
		// //
		// // }
		// // int m = 0;
		// // try{
		// //
		// // m =mc.get(k).getN();
		// //
		// // }catch(NullPointerException ex){
		// //
		// // }
		// // int mu = 0;
		// // try{
		// //
		// // mu =mcu.get(k).getN();
		// //
		// // }catch(NullPointerException ex){
		// //
		// // }
		// //
		// //
		// System.out.println(k+"\t"+amc.get(k).getN()+"\t"+au+"\t"+m+"\t"+mu);
		// //
		// // }

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
				
//				if(chr.equals("6")||chr.equals("8")){
//					str = br.readLine();
//					continue;
//				}
//					
					
//				if(poskey.equals("6-397352")){
//					System.out.println("here"+str);
//				}else{
//					System.out.println(str);
//				}
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
