import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;


public class Main {
	
	public static void benchmark(int ml, int nm, int sl, int sc, int setnum) {
		char[][] seqs = new char[sc][sl];
		for(int j=0; j<sc; j++) {
			for(int i=0; i<sl/4; i++) {
				seqs[j][i] = 'A';
				seqs[j][i + (sl/4)] = 'C';
				seqs[j][i + (2 * sl / 4)] = 'G';
				seqs[j][i + (3 * sl / 4)] = 'T';
			}
		}
		
		for(int j=0; j<sc; j++) {
			for(int i=0; i<sl; i++) {
				for(int k=i; k<sl; k++) {
					int r = (int) (Math.random() * (sl - k));
					char temp = seqs[j][k];
					seqs[j][k] = seqs[j][r];
					seqs[j][r] = temp;
				}
			}
		}
		
		char[] motif = new char[ml];
		for(int i=0; i<ml; i++) {
			switch((int)(Math.random() * 4)) {
				case 0: motif[i] = 'A'; break;
				case 1: motif[i] = 'C'; break;
				case 2: motif[i] = 'G'; break;
				case 3: motif[i] = 'T'; break;
				default: System.out.println("why rand bork");
			}
		}
		for(int i=0; i<nm; i++){
			int r=0;
			do {
				r = (int)(Math.random() * ml);
			} while (motif[r] == '*');
			motif[r] = '*';
		}
		
		char[][] plants = new char[sc][ml];
		for(int j=0; j<sc; j++) {
			for(int i=0; i<ml; i++) {
				plants[j][i] = motif[i];
				if(motif[i] == '*') {
					switch((int)(Math.random() * 4)) {
						case 0: plants[j][i] = 'A'; break;
						case 1: plants[j][i] = 'C'; break;
						case 2: plants[j][i] = 'G'; break;
						case 3: plants[j][i] = 'T'; break;
						default: System.out.println("why rand bork");
					}	
				}
			}
		}
		
		(new File("data/set" + setnum)).mkdir();
		File seqF = new File("data/set"+setnum+"/sequences.fa");
		File sitesF = new File("data/set"+setnum+"/sites.txt");
		File motifF = new File("data/set"+setnum+"/motif.txt");
		File mlenF = new File("data/set"+setnum+"/motiflength.txt");
		
		try {
			seqF.createNewFile();
			sitesF.createNewFile();
			motifF.createNewFile();
			
			PrintWriter seqP = new PrintWriter(new FileWriter(seqF));
			PrintWriter sitesP = new PrintWriter(new FileWriter(sitesF));
			PrintWriter motifP = new PrintWriter(new FileWriter(motifF));
			PrintWriter mlenP = new PrintWriter(new FileWriter(mlenF));
			
			for(int j=0; j<sc; j++) {
				int r = (int)(Math.random() * (sl-ml+1));
				sitesP.println(r);
				for(int i=0; i<ml; i++) {
					seqs[j][r+i] = plants[j][i];
				}
			}
			
			for(int j=0; j<sc; j++) {
				seqP.println(">set" + setnum + "sequence "+j);
				for(int i=0; i<sl; i++) {
					seqP.print(seqs[j][i]);
				}
				seqP.println();
			}
			
			motifP.print("MOTIF" + setnum + "    " + ml + "    ");
			for(int i=0; i<ml; i++) {
				motifP.print(motif[i]);
			}
			
			mlenP.print(ml);
			
			seqP.close();
			sitesP.close();
			motifP.close();
			mlenP.close();
			
			System.out.println("Success on set " + setnum);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
			
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//ml, nm, sl, sc
		for(int i=0; i<10; i++) {
			benchmark(8, 1, 500, 10, i);
			benchmark(8, 0, 500, 10, i+10);
			benchmark(8, 2, 500, 10, i+20);
			benchmark(6, 1, 500, 10, i+30);
			benchmark(7, 1, 500, 10, i+40);
			benchmark(8, 1, 500, 5, i+50);
			benchmark(8, 1, 500, 20, i+60);		
		}

	}

}
