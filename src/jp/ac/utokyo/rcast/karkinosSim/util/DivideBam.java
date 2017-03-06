/*
 * @author Hiroki Ueda
 * ueda[at]genome.rcast.u-tokyo.ac.jp
 * @version 1.0
 * @since 1.0
 */
package jp.ac.utokyo.rcast.karkinosSim.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

public class DivideBam {

	public static void main(String[] arg) {


		try {
			
			//String sam1 = "/GLUSTER_DIST/data/users/ueda/karkinosSim/THCC135/THCC_135T-THCC_135N_normal_pre.bam";
			//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim/THCC135/";
			//String sam1 = "/GLUSTER_DIST/data/users/ueda/karkinosSim/test/C39_normal_genome.bam";
			//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim/test/";
			//exec(sam1,outdir);

			
			//String sam1 = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bam/THCC_12-15T-THCC_12-15N_normal_pre.bam";
			//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bamsim";
			
			
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
				help.printHelp("karkinosSim.jar devideBAM", opts, true);
				return;
			}
			
			String sam1 = cl.getOptionValue("bam");
			String outdir = cl.getOptionValue("out");
					
			exec(sam1,outdir);

//			String sam2 = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bam/THCC_150T-THCC_150N_normal_pre.bam";
//			//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinosSim2/bamsim";
//			exec(sam2,outdir);
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();

		
		optionlist.add(getOption("bam", "bam", true,
				"input bam file", true));

		optionlist.add(getOption("out", "out", true, "output directory", true));

		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	private static void exec(String sam1, String outdir) throws Exception {
		
		SAMFileReader sr1 = getReader(sam1);
		
		SAMFileHeader sfh = sr1.getFileHeader().clone();
		
		File f = new File(sam1);
		String name = f.getName();
		name = name.replaceAll(".bam", "");
		
		SAMFileWriter writer0 = getWriter(sfh, outdir+"/"+name+"_0.bam");
		SAMFileWriter writer1 = getWriter(sfh, outdir+"/"+name+"_1.bam");
		SAMFileWriter writer2 = getWriter(sfh, outdir+"/"+name+"_2.bam");
		
		
		Iterator<SAMRecord> ite = sr1.iterator();
		while(ite.hasNext()){
			
			SAMRecord sam = ite.next();
			String readname = sam.getReadName();
			String ypic = readname.substring(readname.lastIndexOf(":")+1);
			int ypixint = 0;
			try{
				ypixint = Integer.parseInt(ypic);
			}catch(Exception ex){
				
			}
			int mod3 = ypixint%3;
			if(mod3==0){
				writer0.addAlignment(sam);
			}else if(mod3==1){
				writer1.addAlignment(sam);
			}else{
				writer2.addAlignment(sam);
			}

		}
		sr1.close();
		writer0.close();
		writer1.close();
		writer2.close();
	}

	private static SAMFileReader getReader(String in) {

		File INPUT = new File(in);
		SAMFileReader reader = new SAMFileReader(INPUT);
		SAMFileHeader sfh = reader.getFileHeader();
		if (sfh.getAttribute("SO") == null
				|| sfh.getAttribute("SO").equals("sorted")) {
			sfh.setSortOrder(SortOrder.coordinate);
		}
		reader.setValidationStringency(ValidationStringency.SILENT);
		return reader;
	}

	private static SAMFileWriter getWriter(SAMFileHeader sfh2, String OUTPUT)
			throws Exception {


		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		sfh2.setSortOrder(SortOrder.coordinate);
		factory.setCreateIndex(true);
		SAMFileWriter writer = factory.makeSAMOrBAMWriter(sfh2, true, new File(
				OUTPUT));

		return writer;
	}
}
