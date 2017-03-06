package jp.ac.utokyo.rcast.karkinosSim;

import java.io.IOException;
import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class ReadsAdds {

	public static void addsReads(TwoBitGenomeReader tgr, int cn,
			List<SamBean> alleleAReads, List<SamBean> alleleAReadsAdd) {

		//
		if (cn == 2) {

			alleleAReads.addAll(alleleAReadsAdd);

		} else if (cn > 2) {

			// add artifitial reads
			int addnum = cn - 2;
			int serial = 0;
			for (int n = 0; n < addnum; n++) {

				for (SamBean sb : alleleAReadsAdd) {
					serial++;
					// make artifitial amplification
					if (sb.bothExsist()
							&& (sb.sam1.getDuplicateReadFlag() == false)) {
						SamBean sbmod = draft(tgr, sb, serial);
						alleleAReads.add(sbmod);
					}

				}

			}

		}
		//

	}

	public static SamBean draft(TwoBitGenomeReader tgr, SamBean sb, int serial) {

		SamBean sbean = new SamBean();
		String name = sb.sam1.getReadName();
		name = name + ":" + serial;
		// draft 1-30nt
		SAMRecord sam1 = null;
		SAMRecord sam2 = null;
		try {
			sam1 = draft(name, sb.sam1, tgr);
			sam2 = draft(name, sb.sam2, tgr);
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		sbean.add(sam1);
		sbean.add(sam2);
		return sbean;

	}

	public static int drafrrange = 15;

	private static SAMRecord draft(String name, SAMRecord sam1,
			TwoBitGenomeReader tgr) throws CloneNotSupportedException {

		SAMRecord sr = (SAMRecord) sam1.clone();
		List<CigarElement> ce = sr.getCigar().getCigarElements();
		if (ce == null || ce.size() == 0) {
			return sr;
		}
		CigarElement first = ce.get(0);
		boolean noindel = first.getLength() == sr.getReadLength();
		if (noindel) {
			// -15 to 15
			int draftsize = rand(-1 * drafrrange, drafrrange);
			//
			int as = sr.getAlignmentStart();
			int ae = sr.getAlignmentEnd();

			String seq = sr.getReadString();
			try {
				seq = draft(tgr, sr.getReferenceName(), as, ae, draftsize, seq);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			as = as + draftsize;
			ae = ae + draftsize;
			sr.setReadString(seq);
			sr.setAlignmentStart(as);
			// sr.setAlignmentEnd(ae);

		}
		return sr;
	}

	private static String draft(TwoBitGenomeReader tgr, String chr, int as,
			int ae, int draftsize, String seq) throws IOException {

		int asa = as + draftsize;
		int aea = ae + draftsize;
		//
		if (draftsize == 0) {
			return seq;
		}
		int seqlen = seq.length();
		if (draftsize > 0) {

			// left chop and append to end
			String s = seq.substring(draftsize);
			String add = tgr.getGenomicSeq(chr, ae + 1, aea, true);
			return s + add;

		} else {

			// right chop and insert to first
			String s = seq.substring(0, seqlen + draftsize);
			String add = tgr.getGenomicSeq(chr, asa, as - 1, true);
			return add + s;

		}
	}

	private static int rand(int start, int end) {

		//
		int m = Math.random() < 0.5 ? -1 : 1;
		double d = Math.random() * drafrrange * m;
		int n = Math.round((float) d);
		if (n < start)
			n = start;
		if (n > end)
			n = end;
		return n;

	}

	public static SamBean draft(TwoBitGenomeReader tgr, SamBean sb, int serial,
			SNPHolder snpHolder) {

		SamBean sbean = new SamBean();
		String name = sb.sam1.getReadName();
		name = name + ":" + serial;
		// draft 1-30nt
		SAMRecord sam1 = null;
		SAMRecord sam2 = null;
		try {
			sam1 = draft(name, sb.sam1, tgr, snpHolder);
			sam2 = draft(name, sb.sam2, tgr, snpHolder);
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		sbean.add(sam1);
		sbean.add(sam2);
		return sbean;
	}

	private static SAMRecord draft(String name, SAMRecord sam1,
			TwoBitGenomeReader tgr,SNPHolder snpHolder) throws CloneNotSupportedException {

		SAMRecord sr = (SAMRecord) sam1.clone();
		List<CigarElement> ce = sr.getCigar().getCigarElements();
		if (ce == null || ce.size() == 0) {
			return sr;
		}		
		
		CigarElement first = ce.get(0);
		boolean noindel = first.getLength() == sr.getReadLength();
		if (noindel) {
			// -15 to 15
			int draftsize = rand(-1 * drafrrange, drafrrange);
			//
			int as = sr.getAlignmentStart();
			int ae = sr.getAlignmentEnd();

			String seq = sr.getReadString();
			try {
				seq = draft(tgr, sr.getReferenceName(), as, ae, draftsize, seq);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			as = as + draftsize;
			ae = ae + draftsize;
			
			try {
				seq = checkSeq(tgr,seq,sr.getReferenceName(), as, ae,snpHolder);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			sr.setReadString(seq);
			sr.setAlignmentStart(as);
			// sr.setAlignmentEnd(ae);

		}
		return sr;
	}

	private static String checkSeq(TwoBitGenomeReader tgr,String seq, String referenceName, int as,
			int ae, SNPHolder snpHolder) throws IOException {
		
		StringBuffer sb = new StringBuffer();
		int m=0;
		for(int n=as;n<=ae;n++){
			
			char a = tgr.getGenomeNuc(referenceName, n, true);
			char b = seq.charAt(m);
			if(a==b){
				sb.append(a);
			}else{
				
				if(snpHolder.snpmap.containsKey(n)){
					sb.append(b);
				}else{
					sb.append(a);
				}
				
			}
			m++;
		}
		return sb.toString();
		
	}

}
