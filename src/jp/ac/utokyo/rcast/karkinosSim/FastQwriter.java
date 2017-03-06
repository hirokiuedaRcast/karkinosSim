package jp.ac.utokyo.rcast.karkinosSim;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import jp.ac.utokyo.rcast.karkinosSim.util.FastqWriterOs;

public class FastQwriter {

	FastqWriterOs fq1 = null;
	FastqWriterOs fq2 = null;
	long count = 0;
	int serial = 0;
	String outdir;
	String id;
	int unit = 0;

	public FastQwriter(String outdir, String id, int unit) {

		this.outdir = outdir;
		File f = new File(outdir);
		if(!f.exists()){
			f.mkdirs();
		}
		this.id = id;
		this.unit = unit;

	}

	public void close() {

		if (fq1 != null) {
			fq1.close();
		}
		if (fq2 != null) {
			fq2.close();
		}

	}

	private void renewFile() throws IOException {

		close();
		//
		String fqs1 = outdir + "/" + id + "_R1_" + toStr(serial) + ".fastq.gz";
		String fqs2 = outdir + "/" + id + "_R2_" + toStr(serial) + ".fastq.gz";
		OutputStream gzipOutStream1 = new GZIPOutputStream(
				new BufferedOutputStream(new FileOutputStream(fqs1)));
		OutputStream gzipOutStream2 = new GZIPOutputStream(
				new BufferedOutputStream(new FileOutputStream(fqs2)));

		fq1 = new FastqWriterOs(gzipOutStream1);
		fq2 = new FastqWriterOs(gzipOutStream2);

	}

	public void add(SamBean sb) throws IOException {

		if (count % unit == 0) {
			serial++;
			renewFile();
		}
		count++;
		SAMRecord sam1 = sb.sam1;

		String header = sam1.getReadName() + " 1:N:0";
		String readString = sam1.getReadString();
		String bqs = sam1.getBaseQualityString();
		if (sam1.getReadNegativeStrandFlag()) {
			//
			readString = revcon(readString);
			bqs = rev(bqs);
		}
		FastqRecord rec = new FastqRecord(header, readString, "", bqs);
		fq1.write(rec);

		SAMRecord sam2 = sb.sam2;

		header = sam2.getReadName() + " 2:N:0";
		readString = sam2.getReadString();
		bqs = sam2.getBaseQualityString();
		if (sam2.getReadNegativeStrandFlag()) {
			//
			readString = revcon(readString);
			bqs = rev(bqs);
		}
		FastqRecord rec2 = new FastqRecord(header, readString, "", bqs);
		fq2.write(rec2);

	}

	private static String rev(String read) {

		StringBuffer sb = new StringBuffer(read);
		sb.reverse();
		return sb.toString();
	}

	private static String revcon(String read) {

		StringBuffer sb = new StringBuffer();
		for (char c : read.toCharArray()) {
			sb.append(comp(c));
		}
		sb.reverse();
		return sb.toString();
	}

	private static char comp(char c) {

		if (c == 'A') {
			return 'T';
		} else if (c == 'T') {
			return 'A';
		} else if (c == 'C') {
			return 'G';
		} else if (c == 'G') {
			return 'C';
		}
		return 'N';
	}

	private String toStr(int fqserial) {

		String s = String.valueOf(fqserial);
		if (s.length() == 1) {
			return "00" + s;
		} else if (s.length() == 2) {
			return "0" + s;
		} else {
			return s;
		}

	}

}
