package jp.ac.utokyo.rcast.karkinosSim.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedHashSet;
import java.util.Set;

public class FIlterVcf {

	public static void main(String[] arg) {

		 String in
		 ="/GLUSTER_DIST/data/Genomes/karkinos/1000G_biallelic.indels.hg19.vcf";
		 String
		 out="/GLUSTER_DIST/data/users/ueda/karkinosSim/pg/mutect/ref/1000G_biallelic.indels.hg19_2.vcf";

//		String in = "/GLUSTER_DIST/data/users/ueda/karkinosSim/pg/mutect/ref/dbsnp_132_b37.leftAligned.vcf";
//		String out = "/GLUSTER_DIST/data/users/ueda/karkinosSim/pg/mutect/ref/dbsnp_132_b37.leftAlignedmod2.vcf";

		try {

			Set<String> s = new LinkedHashSet<String>();

			BufferedWriter bw = new BufferedWriter(new FileWriter(out));

			BufferedReader br = new BufferedReader(new FileReader(in));
			String str = br.readLine();
			int cnt = 0;
			while (str != null) {
				cnt++;
				if (str.startsWith("#")) {
					bw.write(str + "\n");
					str = br.readLine();
					continue;
				}
				String[] sa = str.split("\t");
				String chr = sa[0];
				chr = chr.replace("chr", "");
				if (chr.equals("MT")) {
					chr = "M";
				}
				String chr2 = "chr" + chr;
				sa[0] = chr2;
				int cntc = 0;
				for (String ss : sa) {

					if (cntc > 0) {
						bw.write("\t");
					}
					bw.write(ss);
					cntc++;
				}
				bw.write("\n");

				s.add(chr);
				str = br.readLine();

			}

			br.close();
			bw.close();

			for (String ss : s) {
				System.out.println(ss);
			}

		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

}
