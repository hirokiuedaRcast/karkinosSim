package jp.ac.utokyo.rcast.karkinosSim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class CNVHolder {

	public CNVHolder(File file) {
		try {
			init(file);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	List<CNVBean> list = new ArrayList<CNVBean>();
	
	public void init(File file) throws IOException{
		
		FileInputStream fis = new FileInputStream(file);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		for (;;) {
			String line = br.readLine();
			if (line == null)break;
			if(line.startsWith("#"))continue;
			list.add(new CNVBean(line));
			
		}	
	}

	public List<CNVBean> getList() {
		return list;
	}

}
