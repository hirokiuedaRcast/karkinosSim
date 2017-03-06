package jp.ac.utokyo.rcast.karkinosSim.util;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintStream;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

public class FastqWriterOs {
	private final File file;
	private final PrintStream writer;
	private OutputStream os;

	public FastqWriterOs(final File file) {
		this.file = file;
		this.writer = new PrintStream(IOUtil.openFileForWriting(file));
	}

	public FastqWriterOs(final OutputStream os) {

		this.os = os;
		this.writer = new PrintStream(os);
		this.file = null;
	}

	public void write(final FastqRecord rec) {
		writer.print(FastqConstants.SEQUENCE_HEADER);
		writer.println(rec.getReadHeader());
		writer.println(rec.getReadString());
		writer.print(FastqConstants.QUALITY_HEADER);
		writer.println(rec.getBaseQualityHeader() == null ? "" : rec
				.getBaseQualityHeader());
		writer.println(rec.getBaseQualityString());
		if (writer.checkError()) {
			throw new PicardException("Error in writing file " + file);
		}
	}

	public void close() {
		writer.close();
		if(os!=null){
			try{
				os.close();
			}catch(Exception ex){
				
			}
		}
		
	}

}