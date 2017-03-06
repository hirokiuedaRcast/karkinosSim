package jp.ac.utokyo.rcast.karkinosSim;

/*
 * Copyright Hiroki Ueda

 *  This program is free software; you can redistribute it and/or modify it under
 *	the terms of the GNU General Public License as published by the Free Software
 *	Foundation; either version 2 of the License, or (at your option) any later
 *	version.
	
 *	This program is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *	FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 *	details.
	
 *	You should have received a copy of the GNU General Public License along with
 *	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 *	Place, Suite 330, Boston, MA 02111-1307 USA
 */


import java.io.IOException;

import jp.ac.utokyo.karkinosSim.candidategen.RomdomAssign;
import jp.ac.utokyo.rcast.karkinosSim.util.AnswerCheck;
import jp.ac.utokyo.rcast.karkinosSim.util.DivideBam;



public class KarkinosSimCmd {

	
	public static final String version = "karkinosSim version 0.1.0 2017/03/06";
	
	public static void main(String[] arg) throws Exception {

		//
		// arg = new String[]{"analysis"};
		if (arg == null || arg.length == 0) {
			printMsg();
		} else {

			String cmd = arg[0];
			String[] arg2 = getArg(arg);
			
			if (cmd.equals("devideBAM")) {

				DivideBam.main(arg2);	
				
			}else if(cmd.equals("assignSNVIndel")){
				
				RomdomAssign.main(arg2);
				
			} else if (cmd.equals("similate")) {

				KarkinosSimMain.main(arg2);	
				
			} else if (cmd.equals("checkAnswer")) {

				AnswerCheck.main(arg2);	
				
			}else{

				printMsg();
			}

		}

	}

	private static String[] getArg(String[] arg) {
		String[] arg2 = new String[arg.length - 1];
		for (int n = 1; n < arg.length; n++) {
			arg2[n - 1] = arg[n];
			// System.out.println(arg[n]);
		}
		return arg2;
	}

	private static void printMsg() {
		System.out.println(version);
		System.out.println("usage: karkinosSim  <command> options");
		System.out.println("possible cmmands are: devideBAM , assignSNVIndel, similate, checkAnswer");
		

	}

}
