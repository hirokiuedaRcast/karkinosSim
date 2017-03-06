package jp.ac.utokyo.rcast.karkinosSim;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class ReadsSimilatorTetro {

	List<SamBean> sambeanList = new ArrayList<SamBean>();

	int smutationExsistFlg = 0;
	CNVBean bean;
	CNVBean beanfocal;
	CloseableIterator<SAMRecord> normarIte1;
	SAMFileReader normalbam2;
	GenomeInterval gi;
	CNVHolder cnv;
	List<SNPInfoBean> artifitialMutation;
	List<SAMRecord> sortedSAMRecords = new ArrayList<SAMRecord>();
	double tc;
	SNPHolder snpHolder;
	TwoBitGenomeReader tgr;

	public ReadsSimilatorTetro(TwoBitGenomeReader tgr, SNPHolder snpHolder,
			GenomeInterval gi, CNVHolder cnv, double tc,
			List<SNPInfoBean> artifitialMutation, int smutationExsistFlg,
			CNVBean bean, CNVBean cnvbeanfocal,
			CloseableIterator<SAMRecord> normarIte1, SAMFileReader normalbam2) {

		//
		this.gi = gi;
		this.cnv = cnv;
		this.artifitialMutation = overlap(artifitialMutation);
		this.smutationExsistFlg = smutationExsistFlg;
		this.bean = bean;
		this.normarIte1 = normarIte1;
		this.normalbam2 = normalbam2;
		this.tc = tc;
		this.snpHolder = snpHolder;
		beanfocal = cnvbeanfocal;
		this.tgr = tgr;
		this.sambeanList = new ArrayList<SamBean>();
	}

	private List<SNPInfoBean> overlap(List<SNPInfoBean> artifitialMutation2) {

		List<SNPInfoBean> list = new ArrayList<SNPInfoBean>();
		//
		for (SNPInfoBean sb : artifitialMutation2) {

			//
			if (sb.chr.equals(gi.chr)) {

				//
				if ((gi.start <= sb.pos) && (gi.end >= sb.pos)) {
					list.add(sb);
				}

			}

		}
		return list;
	}

	public void modifyReads() throws IOException {

		// divide into nomodyfy reads for normal contamination,allele A,allele B
		List<SamBean> keepReads = new ArrayList<SamBean>();

		List<SamBean> alleleA1Reads = new ArrayList<SamBean>();
		List<SamBean> alleleA2Reads = new ArrayList<SamBean>();
		List<SamBean> alleleA3Reads = new ArrayList<SamBean>();
		List<SamBean> alleleA4Reads = new ArrayList<SamBean>();

		List<SamBean> alleleB1Reads = new ArrayList<SamBean>();
		List<SamBean> alleleB2Reads = new ArrayList<SamBean>();
		List<SamBean> alleleB3Reads = new ArrayList<SamBean>();
		List<SamBean> alleleB4Reads = new ArrayList<SamBean>();

		CloseableIterator<SAMRecord> normarIte2 = normalbam2.query(gi.getChr(),
				gi.getStart(), gi.getEnd(), false);

		devideDataInto4(normarIte1, normarIte2, keepReads, alleleA1Reads,
				alleleA2Reads, alleleA3Reads, alleleA4Reads, alleleB1Reads,
				alleleB2Reads, alleleB3Reads, alleleB4Reads);

		//
		List<SamBean> alleleAReads = new ArrayList<SamBean>();
		List<SamBean> alleleBReads = new ArrayList<SamBean>();

		//
		if (bean.cnA >= 1) {
			alleleAReads.addAll(alleleA1Reads);
		}
		if (bean.cnA >= 2) {
			alleleAReads.addAll(alleleA2Reads);
		}
		if (bean.cnA >= 3) {
			alleleAReads.addAll(alleleA3Reads);
		}
		if (bean.cnA >= 4) {
			alleleAReads.addAll(alleleA4Reads);
		}

		if (bean.cnA > 4) {

			addArtifitial(alleleAReads, alleleA1Reads, bean.cnA);

		}

		if (bean.cnB >= 1) {
			alleleBReads.addAll(alleleB1Reads);
		}
		if (bean.cnB >= 2) {
			alleleBReads.addAll(alleleB2Reads);
		}
		if (bean.cnB >= 3) {
			alleleBReads.addAll(alleleB3Reads);
		}
		if (bean.cnB >= 4) {
			alleleBReads.addAll(alleleB4Reads);
		}

		//
		induceMutation(tgr, bean.cnA, alleleAReads);

		// combine reads
		combine(keepReads, alleleAReads, alleleBReads);
		// debug
		// sort by pos
		sortBYpos(); // no

	}

	private void addArtifitial(List<SamBean> alleleAReads,
			List<SamBean> alleleA1Reads, int cnA) {

		int add = (bean.cnA - 4);
		int serial = 0;

		while (add > 0) {

			for (SamBean sb : alleleA1Reads) {
				serial++;
				if (sb.bothExsist()
						&& (sb.sam1.getDuplicateReadFlag() == false)) {
					SamBean sbmod = ReadsAdds.draft(tgr, sb, serial, snpHolder);
					alleleAReads.add(sbmod);
				}
			}
			add--;
		}

	}

	private void devideDataInto4(CloseableIterator<SAMRecord> normarIte12,
			CloseableIterator<SAMRecord> normarIte2, List<SamBean> keepReads,
			List<SamBean> alleleA1Reads, List<SamBean> alleleA2Reads,
			List<SamBean> alleleA3Reads, List<SamBean> alleleA4Reads,
			List<SamBean> alleleB1Reads, List<SamBean> alleleB2Reads,
			List<SamBean> alleleB3Reads, List<SamBean> alleleB4Reads) {

		//
		LinkedHashMap<String, SamBean> readsmap = new LinkedHashMap<String, SamBean>();
		while (normarIte1.hasNext()) {

			SAMRecord sr = normarIte1.next();
			add(readsmap, sr);
		}
		normarIte1.close();
		Set<Entry<String, SamBean>> set = readsmap.entrySet();

		for (Entry<String, SamBean> et : set) {

			//
			// String rn = et.getKey();

			SamBean sb = et.getValue();
			boolean keepAsNormalContamination = keepAsNormalContamination(tc);
			if (keepAsNormalContamination) {

				keepReads.add(sb);

			} else {

				boolean aallele = assignAorB(sb);
				boolean rand = Math.random() < 0.5;
				if (aallele) {

					if (rand) {
						alleleA1Reads.add(sb);
					} else {
						alleleA2Reads.add(sb);
					}

				} else {

					if (rand) {
						alleleB1Reads.add(sb);
					} else {
						alleleB2Reads.add(sb);
					}
				}
			}
		}
		// balance
		// 2n is assumed, so if A and B are out of balance ,adjusted
		// balance(alleleAReadslocal,alleleBReadslocal);
		readsmap = new LinkedHashMap<String, SamBean>();
		while (normarIte2.hasNext()) {

			SAMRecord sr = normarIte2.next();
			add(readsmap, sr);
		}
		normarIte2.close();
		set = readsmap.entrySet();

		for (Entry<String, SamBean> et : set) {

			SamBean sb = et.getValue();
			boolean keepAsNormalContamination = keepAsNormalContamination(tc);
			if (keepAsNormalContamination) {

				// keepReads.add(sb);

			} else {

				boolean aallele = assignAorB(sb);
				boolean rand = Math.random() < 0.5;
				if (aallele) {

					if (rand) {
						alleleA3Reads.add(sb);
					} else {
						alleleA4Reads.add(sb);
					}

				} else {

					if (rand) {
						alleleB3Reads.add(sb);
					} else {
						alleleB4Reads.add(sb);
					}
				}

			}
		}

	}

	private char getNuc(SAMRecord sam, int pos) {

		if (sam == null)
			return 'N';
		IndelInfo ii = new IndelInfo();
		return AssinUtils.getNuc(pos, sam, ii);
	}

	private void induceMutation(TwoBitGenomeReader tgr, int cnA,
			List<SamBean> alleleAReads) throws IOException {

		//
		for (SNPInfoBean sib : artifitialMutation) {
			//
			int pos = sib.pos;
			if (Math.abs(pos - 93726528) < 100) {
				System.out.println("stop here");
			}
			for (SamBean sb : alleleAReads) {

				boolean domodify = doModfy(cnA, pos);
				if (domodify) {
					if (induceErrorCoffef(0.99)) {
						domodify = false;
					}
					// subpolulation ratio 0.2
					if (sib.modeSnp == 3) {
						if (induceErrorCoffef(0.5)) {
							domodify = false;
						}
					}

				}
				if (domodify) {
					boolean c1 = contain(pos, sb.sam1);
					if (c1) {
						//
						modify(tgr, sib, sb.sam1);
					}

					boolean c2 = contain(pos, sb.sam2);
					if (c2) {
						//
						modify(tgr, sib, sb.sam2);
					}
				}

			}

		}

	}

	private boolean induceErrorCoffef(double d) {

		return Math.random() > d;
	}

	private void modify(TwoBitGenomeReader tgr2, SNPInfoBean sib, SAMRecord sam)
			throws IOException {

		IndelInfo ii = new IndelInfo();
		int pos = sib.pos;
		int mutationidx = KarkinosSimUtils.getCharIdx(pos, sam, ii);
		if (mutationidx < 0)
			return;
		String ref = sib.ref;
		String alt = sib.alt;

		if (sib.isIndel()) {

			// do not performe cigar collection since
			// reads are mapped by aligner again anyway
			int reflen = ref.length() - 1;
			int altlen = alt.length() - 1;
			String rs = sam.getReadString();
			if (reflen > altlen) {

				String add = tgr.getGenomicSeq(sam.getReferenceName(),
						sam.getAlignmentEnd() + 1, sam.getAlignmentEnd()
								+ reflen, true);
				// deletion
				StringBuffer sb = new StringBuffer();
				for (int n = 0; n < rs.length(); n++) {

					//
					if (n >= mutationidx && n < (mutationidx + reflen)) {

					} else {
						sb.append(rs.charAt(n));
					}

				}
				sb.append(add);
				int seqlen = sam.getReadLength();
				String seqo = substring(sb.toString(), 0, seqlen);
				sam.setReadString(seqo);

			} else {

				// insersion
				int seqlen = sam.getReadLength();
				StringBuffer sb = new StringBuffer();
				for (int n = 0; n < rs.length(); n++) {

					//
					if (n == mutationidx) {
						sb.append(alt.substring(1));
					} else {
						sb.append(rs.charAt(n));
					}

				}
				String seqo = substring(sb.toString(), 0, seqlen);
				sam.setReadString(seqo);

			}

		} else {

			//

			byte b = sam.getBaseQualities()[mutationidx];
			boolean error = error(b);
			if (error) {
				alt = randbase();
			}
			byte[] seq = sam.getReadBases();
			seq[mutationidx] = alt.getBytes()[0];
			sam.setReadBases(seq);

		}

	}

	private String substring(String str, int s, int e) {

		return KarkinosSimUtils.substring(str, s, e);
	}

	private String randbase() {

		double rand = Math.random() * 4;
		if ((0 <= rand) && (1 > rand)) {
			return "A";
		} else if ((1 <= rand) && (2 > rand)) {
			return "T";
		} else if ((2 <= rand) && (3 > rand)) {
			return "G";
		} else {
			return "C";
		}
	}

	private boolean error(byte b) {

		double qual0 = (int) b & 0xFF;
		qual0 = qual0 * 0.1;
		double pNomatch = (1 / Math.pow(10, qual0));
		return Math.random() < pNomatch;
	}

	private boolean doModfy(int cnA, int pos) {

		if (cnA == 2) {

			if (pos % 2 == 0) {
				return true;

			} else {

				return Math.random() < 0.5;
			}

		}
		if (cnA == 3) {

			if (pos % 3 == 0) {
				return true;
			} else if (pos % 3 == 1) {
				return Math.random() < 0.33;
			} else {
				return Math.random() < 0.66;
			}
		}
		return true;
	}

	private boolean contain(int pos, SAMRecord sr) {

		if (sr != null) {
			if (pos >= sr.getAlignmentStart() && pos <= sr.getAlignmentEnd()) {
				return true;
			}
		}
		return false;
	}

	private void sortBYpos() {

		TreeMap<Integer, SAMRecord> tm = new TreeMap<Integer, SAMRecord>();
		for (SamBean sb : sambeanList) {

			if (sb.sam1 != null) {
				tm.put(sb.sam1.getAlignmentStart(), sb.sam1);
			}
			if (sb.sam2 != null) {
				tm.put(sb.sam2.getAlignmentStart(), sb.sam2);
			}

		}

		//
		Set<Entry<Integer, SAMRecord>> set = tm.entrySet();
		for (Entry<Integer, SAMRecord> et : set) {
			sortedSAMRecords.add(et.getValue());
		}

		// for (SamBean sb : sambeanList) {
		//
		// if (sb.sam1 != null) {
		// sortedSAMRecords.add(sb.sam1);
		// }
		// if (sb.sam2 != null) {
		// sortedSAMRecords.add(sb.sam2);
		// }
		//
		// }

	}

	private void combine(List<SamBean> keepReads, List<SamBean> alleleAReads,
			List<SamBean> alleleBReads) {

		sambeanList.addAll(keepReads);
		sambeanList.addAll(alleleAReads);
		sambeanList.addAll(alleleBReads);

	}

	private void devideData(CloseableIterator<SAMRecord> normarIte1,
			List<SamBean> keepReads, List<SamBean> alleleAReads,
			List<SamBean> alleleBReads, boolean addnormal, boolean addA,
			boolean addB) {

		LinkedHashMap<String, SamBean> readsmap = new LinkedHashMap<String, SamBean>();
		while (normarIte1.hasNext()) {

			SAMRecord sr = normarIte1.next();
			add(readsmap, sr);
		}
		normarIte1.close();
		Set<Entry<String, SamBean>> set = readsmap.entrySet();

		List<SamBean> alleleAReadslocal = new ArrayList<SamBean>();
		List<SamBean> alleleBReadslocal = new ArrayList<SamBean>();

		for (Entry<String, SamBean> et : set) {

			//
			// String rn = et.getKey();

			SamBean sb = et.getValue();
			boolean keepAsNormalContamination = keepAsNormalContamination(tc);
			if (keepAsNormalContamination) {

				if (addnormal) {
					keepReads.add(sb);
				}

			} else {

				boolean aallele = assignAorB(sb);
				if (aallele) {
					alleleAReadslocal.add(sb);
				} else {
					alleleBReadslocal.add(sb);
				}
			}
		}
		// balance
		// 2n is assumed, so if A and B are out of balance ,adjusted
		// balance(alleleAReadslocal,alleleBReadslocal);

		// add
		if (addA) {
			alleleAReads.addAll(alleleAReadslocal);
		}
		if (addB) {
			alleleBReads.addAll(alleleBReadslocal);
		}
		readsmap = null;
	}

	private void balance(List<SamBean> alleleAReadslocal,
			List<SamBean> alleleBReadslocal) {
		//
		int alAnum = alleleAReadslocal.size() + 1;
		int alBnum = alleleBReadslocal.size() + 1;
		System.out.println(alAnum + "\t" + alBnum);
		int total = alAnum + alBnum + 1;
		double r = (double) alAnum / (double) total;
		if (r > 0.6 || r < 0.4) {
			int diff = Math.abs(alAnum - alBnum);
			double diffr = Math.abs(0.5 - r);
			diffr = diffr * 0.8;
			int move = (int) (diffr * total);
			if (r < 0.4) {
				moveData(alleleBReadslocal, alleleAReadslocal, move);
			} else {
				moveData(alleleAReadslocal, alleleBReadslocal, move);
			}

		}

	}

	private void moveData(List<SamBean> from, List<SamBean> to, int move) {

		int cnt = 0;
		if (move <= 0)
			return;
		for (int n = from.size() - 1; n >= 0; n--) {

			SamBean sb = from.get(n);
			to.add(sb);
			from.remove(n);
			cnt++;
			if (cnt >= move) {
				break;
			}

		}

	}

	private void setABData(CloseableIterator<SAMRecord> normarIte1,
			List<SamBean> alleleAReads, List<SamBean> alleleBReads) {

		LinkedHashMap<String, SamBean> readsmap = new LinkedHashMap<String, SamBean>();
		while (normarIte1.hasNext()) {

			SAMRecord sr = normarIte1.next();
			add(readsmap, sr);
		}
		normarIte1.close();
		Set<Entry<String, SamBean>> set = readsmap.entrySet();
		for (Entry<String, SamBean> et : set) {

			//
			SamBean sb = et.getValue();
			String rn = et.getKey();

			boolean keepAsNormalContamination = keepAsNormalContamination(tc);
			if (keepAsNormalContamination) {

			} else {

				boolean aallele = assignAorB(sb);
				if (aallele) {
					alleleAReads.add(sb);
				} else {
					alleleBReads.add(sb);
				}
			}
		}
		readsmap = null;
	}

	private boolean assignAorB(SamBean sb) {
		return AssinUtils.assignAorB(snpHolder.afmap, sb);
	}

	private boolean keepAsNormalContamination(double tc2) {

		double rand = Math.random();
		boolean tumorread = (rand <= tc);
		return !tumorread;
	}

	public List<SAMRecord> getSortedSAMRecords() {
		return sortedSAMRecords;
	}

	public List<SamBean> getSambeanList() {
		return sambeanList;
	}

	private static void add(LinkedHashMap<String, SamBean> map, SAMRecord sr) {

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

	}

}
