package jp.ac.utokyo.rcast.karkinosSim;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;



public class AssinUtils {

	public static boolean assignAorB(Map<Integer, AlleleFragment> afmap,
			SamBean sb) {

		//
		TreeMap<Integer, AlleleFragment> afTmap = (TreeMap<Integer, AlleleFragment>) afmap;

		// if(afTmap.containsKey(899000)){
		// System.out.println("true");
		// }
		if (sb.getSam1() != null) {
			if (sb.getSam1().getReadName().contains("50510")
					||sb.getSam1().getReadName().contains("86535")
					||sb.getSam1().getReadName().contains("100254")
					||sb.getSam1().getReadName().contains("90627")) {
				System.out.println("here");
			}
		}
		//
		int margin = 10000;
		int min = sb.getMin() - margin;
		int max = sb.getMax();
		Map<Integer, AlleleFragment> tailmap = afTmap.tailMap(min);

		Iterator<Integer> iteaf = tailmap.keySet().iterator();
		AlleleFragment afg = null;
		while (iteaf.hasNext()) {

			int key = iteaf.next();
			afg = afmap.get(key);
			if (afg.start <= sb.getMax() && sb.getMin() <= afg.end) {
				break;
			}
			if (sb.getMax() < afg.start) {
				afg = null;
				break;
			}

		}
		if (afg == null) {
			// no clue random assgin
			return getRand();

		} else {
			// look up allele list
			Iterator<Integer> ite = afg.getAllelemap().keySet().iterator();
			while (ite.hasNext()) {
				int pos = ite.next();

				AllelleBase ab = afg.getAllelemap().get(pos);

				boolean poscovered1 = cov(pos, sb.sam1);
				boolean poscovered2 = cov(pos, sb.sam2);
				boolean poscovered = poscovered1 || poscovered2;

				if (poscovered) {

					SAMRecord sam = poscovered1 ? sb.sam1 : sb.sam2;

					String allelA = ab.alleleA;
					String allelB = ab.alleleB;
					IndelInfo ii = new IndelInfo();
					String nuc = getNuc(pos, sam, ii) + "";

					if (ab.indel()) {

						if (ii.indel) {
							return allelA.equals("Indel");
						} else {
							return allelA.equals("ref");
						}
					} else {
						if(pos == 136561557){
							System.out.println(pos+"\t"+ab.alleleA+"-"+ab.alleleB+"\t"+nuc+"\t"+nuc.equals(allelA));
						}
						if (nuc.equals(allelA)) {
							return true;
						} else if (nuc.equals(allelB)) {
							return false;
						}

					}

				}

			}

			// cannot assign
			return getRand();
		}

	}

	private static boolean cov(int pos, SAMRecord sam) {

		if (sam == null)
			return false;
		return sam.getAlignmentStart() <= pos && pos <= sam.getAlignmentEnd();

	}

	private static boolean getRand() {
		// TODO Auto-generated method stub
		return Math.random() < 0.5;
	}

	public static char getNuc(int pos, SAMRecord sam, IndelInfo ii) {

		if (sam == null)
			return 'N';
		int mutationidx = getCharIdx(pos, sam, ii);
		int readslen = sam.getReadLength();
		if ((mutationidx >= 0) && (mutationidx < readslen)) {
			char ch = (char) sam.getReadBases()[mutationidx];
			return ch;
			//
		}
		return 'N';
	}

	private static int getCharIdx(int pos, SAMRecord sam, IndelInfo indelinfo) {

		int start = sam.getAlignmentStart();
		int relpos = pos - start;

		if (relpos == 0) {
			List<CigarElement> list0 = sam.getCigar().getCigarElements();
			if (list0 != null && list0.size() >= 0) {
				CigarElement ce = list0.get(0);
				if (ce.getOperator() == (CigarOperator.SOFT_CLIP)) {
					return ce.getLength();
				}
			} else {
				return 0;
			}

		}

		int readidx = 0;
		int refidx = 0;

		List<CigarElement> list = sam.getCigar().getCigarElements();
		//
		if (list == null || list.size() == 0) {
			// Cigar does not exist which should not happen
			// assgin
			list = new ArrayList<CigarElement>();
			list.add(new CigarElement(sam.getReadLength(), CigarOperator.M));
		}
		int l = 0;
		for (CigarElement ce : list) {

			int len = ce.getLength();
			if (len == sam.getReadLength()) {
				return relpos;
			}

			if (ce.getOperator().consumesReferenceBases()
					&& ce.getOperator().consumesReadBases()) {

				if (relpos <= refidx + len) {

					int readidxr = readidx + (relpos - refidx);
					// check if insersion exsist in next cigar
					if (relpos == refidx + len) {
						if (l + 1 < list.size()) {
							CigarElement nextcigar = list.get(l + 1);
							if (nextcigar.getOperator() == CigarOperator.INSERTION) {
								indelinfo.indel = true;
								indelinfo.length = nextcigar.getLength();
								indelinfo.insersion = substring(
										sam.getReadString(), readidxr, readidxr
												+ indelinfo.length);
								indelinfo.refpos = refidx + len;
							} else if (nextcigar.getOperator() == CigarOperator.DELETION) {
								indelinfo.indel = true;
								indelinfo.length = nextcigar.getLength();
								indelinfo.refpos = refidx + len;
							}
						}
					}
					return readidxr;
				}
				refidx += len;
				readidx += len;

			} else if (ce.getOperator().consumesReferenceBases()) {

				if (ce.getOperator() == CigarOperator.DELETION) {
					if (relpos == refidx + len) {
						indelinfo.indel = true;
						indelinfo.length = len;
						return -1;
					} else if (relpos < refidx + len) {
						return -1;
					}
				}
				refidx += len;
			} else {

				readidx += len;

			}
			l++;
		}
		return readidx;

	}

	private static String substring(String str, int s, int e) {

		if (e >= str.length()) {
			return str.substring(s);
		} else {
			return str.substring(s, e);
		}

	}

}
