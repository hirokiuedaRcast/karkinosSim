package jp.ac.utokyo.rcast.karkinosSim;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class SNPHolder {

	Map<Integer, SNPInfoBean> snpmap = new TreeMap<Integer, SNPInfoBean>();

	public void add(String chr, int pos, float af, String id, String ref,
			String alt, int modeSnp) {

		if (!snpmap.containsKey(pos)) {
			snpmap.put(pos,
					new SNPInfoBean(chr, pos, af, id, ref, alt, modeSnp));
		}

	}

	Map<Integer, AlleleFragment> afmap = new TreeMap<Integer, AlleleFragment>();

	public int checkLink(int librarysizemax) {

		//
		Iterator<Integer> ite = snpmap.keySet().iterator();
		int linkcnt = 0;
		SNPInfoBean b4 = null;
		int lastpos = 0;
		AlleleFragment lastflg = null;
		while (ite.hasNext()) {

			SNPInfoBean snp = snpmap.get(ite.next());
			
//			if(snp.pos==899000){
//				System.out.println("here");
//			}
			
			if (overlap(b4, snp, librarysizemax)) {
				linkcnt++;
				b4.chainPos = snp.pos;
				snp.chainPos = b4.pos;

				AlleleFragment flgment = null;
				if (lastpos == b4.pos) {
					flgment = lastflg;
					flgment.add(snp);
				} else {
					flgment = new AlleleFragment();
					flgment.add(b4);
					flgment.add(snp);
					afmap.put(b4.pos, flgment);
				}
				lastpos = snp.pos;
				lastflg = flgment;
			} else {
				if (b4 != null) {
					AlleleFragment flgment = new AlleleFragment();
					flgment.add(b4);
					afmap.put(b4.pos, flgment);
				}
			}
			b4 = snp;
		}
		return linkcnt;

	}

	public Map<Integer, SNPInfoBean> getSnpmap() {
		return snpmap;
	}

	public Map<Integer, AlleleFragment> getAfmap() {
		return afmap;
	}

	private boolean overlap(SNPInfoBean b4, SNPInfoBean snp, int librarysizemax) {

		if (b4 == null || snp == null) {
			return false;
		}
		int dist = Math.abs(b4.pos - snp.pos);

		return dist < librarysizemax;
	}

}
