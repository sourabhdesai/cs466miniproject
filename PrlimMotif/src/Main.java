import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;


public class Main {

    public static void benchmark(int ml, int nm, int sl, int sc, int setnum) {

        //Create a base strings with equal A, C, G, T
        char[][] seqs = new char[sc][sl];
        for(int j=0; j<sc; j++) {
            for(int i=0; i<sl/4; i++) {
                seqs[j][i] = 'A';
                seqs[j][i + (sl/4)] = 'C';
                seqs[j][i + (2 * sl / 4)] = 'G';
                seqs[j][i + (3 * sl / 4)] = 'T';
            }
        }

        //Random shuffles
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

        //Create random motif string
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

        //Set random positions as variable
        //Exclude the ends because this essentially is a just a smaller length motif
        //Becomes unfair for the algo to guess with end the * is at.
        for(int i=0; i<nm; i++){
            int r=0;
            do {
                r = (int)(Math.random() * (ml-2));
            } while (motif[r+1] == '*');
            motif[r+1] = '*';
        }

        //Create motifs for each sequence (filling in the * with a random letter)
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
            mlenF.createNewFile();

            PrintWriter seqP = new PrintWriter(new FileWriter(seqF));
            PrintWriter sitesP = new PrintWriter(new FileWriter(sitesF));
            PrintWriter motifP = new PrintWriter(new FileWriter(motifF));
            PrintWriter mlenP = new PrintWriter(new FileWriter(mlenF));

            //Copy over motifs to sequences at random position
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

            //System.out.println("Success on set " + setnum);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void findMotif(int setnum) {
        File seqF = new File("data/set"+setnum+"/sequences.fa");
        File mlenF = new File("data/set"+setnum+"/motiflength.txt");
        File timeF = new File("data/set"+setnum+"/runtime.txt");


        try {
            BufferedReader seqR =  new BufferedReader(new FileReader(seqF));
            BufferedReader mlenR =  new BufferedReader(new FileReader(mlenF));

            //Start runtime timer
            long startTime = System.currentTimeMillis();

            //Read in motif length
            int mlen = Integer.parseInt(mlenR.readLine());
            int lines=0;

            //Read in sequences
            while(seqR.readLine() != null) {
                lines++;
            }
            String[] seqs = new String[lines/2];
            seqR.close();
            seqR = new BufferedReader(new FileReader(seqF));
            int i=0;
            while(seqR.readLine() != null) {
                seqs[i]=seqR.readLine();
                i++;
            }

            //Begin algo
            int slen = mlen;
            int score;
            //Hack to speed up run time. Ignores and lmers that don't have atleast 5 matches
            int basescore = 5;
            String sseq;
            ArrayList<lmer> lmers = new ArrayList<lmer>();
            //Compare each position at the first 2 sequences for similarities
            for(i=0; i<500-slen+1; i++) {
                for(int j=0; j<500-slen+1; j++) {
                    score = 0;
                    sseq = "";
                    //Score each combination of positions
                    for(int k=0; k<slen; k++) {
                        if(seqs[0].charAt(i+k) == seqs[1].charAt(j+k)) {
                            score ++;
                            sseq += seqs[0].charAt(i+k);
                        }
                        else {
                            sseq += "*";
                        }
                    }
                    //Save lmers that score above basescore
                    if(score >= basescore) {
                        lmer l = new lmer(sseq);
                        l.score = score;
                        l.loc = new int[seqs.length];
                        l.loc[0] = i;
                        l.loc[1] = j;
                        lmers.add(l);
                    }
                }
            }

            ArrayList<lmer> toRemove;
            //Compare saved lmers with other sequences
            for(int k=2; k<seqs.length; k++) {
                for(i=0; i<500-slen+1; i++) {
                    for(lmer l : lmers) {
                        score = 0;
                        //Score as before
                        for(int j=0; j<slen; j++) {
                            if(seqs[k].charAt(i+j) == l.seq.charAt(j)) {
                                score ++;
                            }
                        }
                        if(score > l.c_score) {
                            l.c_score = score;
                            l.loc[k] = i;
                        }
                    }
                }
                toRemove = new ArrayList<lmer>();
                for(lmer l : lmers) {
                    //Save score of lmers that passed basescore
                    if(l.c_score >= basescore) {
                        l.score += l.c_score;
                        l.c_score = 0;
                        for(int j=0; j<slen; j++) {
                            //Set differences as a *
                            if(seqs[k].charAt(l.loc[k]+j) != l.seq.charAt(j)) {
                                l.seq.setCharAt(j, '*');
                            }
                        }
                    }
                    //Toss those that didn't
                    else {
                        toRemove.add(l);
                    }
                }
                for(lmer l : toRemove) {
                    lmers.remove(l);
                }
            }

            //Sort by score
            Collections.sort(lmers);
            StringBuilder motif;
            if(lmers.size() > 0) {
                motif = new StringBuilder(lmers.get(0).seq);
            }
            else {
                motif = new StringBuilder("");
            }
            //System.out.println(setnum + ": " + motif);

            //End timer
            long endTime   = System.currentTimeMillis();
            long totalTime = endTime - startTime;

            timeF.createNewFile();
            PrintWriter timeP = new PrintWriter(new FileWriter(timeF));
            timeP.println(totalTime);

            seqR.close();
            mlenR.close();
            timeP.close();

            if(lmers.size() > 0) {
                saveOutput(seqs, motif.toString(), lmers.get(0).loc, setnum);
            }
            else {
                saveOutput(null, null, null, setnum);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void saveOutput(String[] seqs, String pmotif, int[] sites, int setnum) {
        File psitesF = new File("data/set"+setnum+"/predictedsites.txt");
        File pmotifF = new File("data/set"+setnum+"/predictedmotif.txt");

        try {
            psitesF.createNewFile();
            pmotifF.createNewFile();

            PrintWriter psitesP = new PrintWriter(new FileWriter(psitesF));
            PrintWriter pmotifP = new PrintWriter(new FileWriter(pmotifF));

            if(sites == null) {
                psitesP.close();
                pmotifP.close();
                return;
            }

            for(int i=0; i<sites.length; i++) {
                psitesP.println(sites[i]);
            }

            //Create PWM
            pmotifP.println(">PMOTIF" + setnum + "    " + pmotif.length());
            for(int i=0; i<pmotif.length(); i++) {
                int as=0;
                int cs=0;
                int gs=0;
                int ts=0;
                for(int j=0; j<seqs.length; j++) {
                    switch(seqs[j].charAt(i+sites[j])) {
                        case 'A': as++; break;
                        case 'C': cs++; break;
                        case 'G': gs++; break;
                        case 'T': ts++; break;
                    }

                }
                pmotifP.println(as + "    " + cs + "    " + gs + "    " + ts + "    ");
            }
            pmotifP.println("<");
            //System.out.println(setnum + " saved");


            psitesP.close();
            pmotifP.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    public static void evaluate(int setnum) {
        File seqF = new File("data/set"+setnum+"/sequences.fa");
        File sitesF = new File("data/set"+setnum+"/sites.txt");
        File motifF = new File("data/set"+setnum+"/motif.txt");
        File mlenF = new File("data/set"+setnum+"/motiflength.txt");
        File psitesF = new File("data/set"+setnum+"/predictedsites.txt");
        File pmotifF = new File("data/set"+setnum+"/predictedmotif.txt");
        File timeF = new File("data/set"+setnum+"/runtime.txt");

        try {
            BufferedReader seqR =  new BufferedReader(new FileReader(seqF));
            BufferedReader sitesR =  new BufferedReader(new FileReader(sitesF));
            BufferedReader mlenR =  new BufferedReader(new FileReader(mlenF));
            BufferedReader psitesR =  new BufferedReader(new FileReader(psitesF));
            BufferedReader pmotifR =  new BufferedReader(new FileReader(pmotifF));
            BufferedReader timeR =  new BufferedReader(new FileReader(timeF));

            if(pmotifR.readLine() == null) {
                return;
            }

            int lines=0;
            while(seqR.readLine() != null) {
                lines++;
            }
            String[] seqs = new String[lines/2];
            seqR.close();
            seqR = new BufferedReader(new FileReader(seqF));
            int x=0;
            while(seqR.readLine() != null) {
                seqs[x]=seqR.readLine();
                x++;
            }

            int mlen = Integer.parseInt(mlenR.readLine());
            int[][] apwm = new int[mlen][4];
            int[][] ppwm = new int[mlen][4];

            String[] row;
            for(int i=0; i<mlen; i++) {
                row = pmotifR.readLine().split(" +");
                for(int j=0; j<4; j++) {
                    ppwm[i][j] = Integer.parseInt(row[j]);
                }
            }

            //Create PWM
            int[] asites = new int[seqs.length];
            for(int i=0; i<seqs.length; i++) {
                asites[i] = Integer.parseInt(sitesR.readLine());
            }
            for(int i=0; i<mlen; i++) {
                apwm[i][0]=0;
                apwm[i][1]=0;
                apwm[i][2]=0;
                apwm[i][3]=0;
                for(int j=0; j<seqs.length; j++) {
                    switch(seqs[j].charAt(i+asites[j])) {
                        case 'A': apwm[i][0]++; break;
                        case 'C': apwm[i][1]++; break;
                        case 'G': apwm[i][2]++; break;
                        case 'T': apwm[i][3]++; break;
                    }
                }
            }

            int[] psites = new int[seqs.length];
            for(int i=0; i<seqs.length; i++) {
                psites[i] = Integer.parseInt(psitesR.readLine());
            }

            //Calculate Relative Entropy
            double re = 0;
            for(int k=0; k<mlen; k++) {
                double rek = 0;
                for(int i=0; i<4; i++) {

                    double pk = (ppwm[k][i] + 0.0) / 10;
                    double qk = (apwm[k][i] + 0.0) / 10;
                    if(pk == 0) {
                        pk = 0.0000000000000001;
                    }
                    if(qk == 0) {
                        qk = 0.0000000000000001;
                    }
                    rek += pk * (Math.log(pk/qk) / Math.log(2.0));
                }
                re += rek;
            }

            int overlap_sites = 0;
            for(int i=0; i<seqs.length; i++) {
                if(asites[i] + mlen > psites[i] && psites[i] + mlen > asites[i]) {
                    overlap_sites ++;
                }
            }
            /*
            System.out.println("Evaluation for Set " + setnum);
            System.out.println("Relative Entropy: " + re);
            System.out.println("Number of overlap regions: " + overlap_sites);
            System.out.println("Runtime in miliseconds: " + Long.parseLong(timeR.readLine()));
            */

            //JSON Print statements
            System.out.println("\t{");
            System.out.println("\t\t\"set\": " + setnum + ',');
            System.out.println("\t\t\"relativeEntropy\": " + re+ ',');
            System.out.println("\t\t\"overlapSites\": " + overlap_sites + ',');
            System.out.println("\t\t\"runtime\": " +  Long.parseLong(timeR.readLine()) );
            System.out.println("\t}" + ( setnum == 69 ? "" : "," ) );


            seqR.close();
            sitesR.close();
            timeR.close();
            mlenR.close();
            psitesR.close();
            pmotifR.close();



        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (NumberFormatException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * @param args
     */
    public static void main(String[] args) {
		for(int i=0; i<10; i++) {
			benchmark(8, 1, 500, 10, i);
			benchmark(8, 0, 500, 10, i+10);
			benchmark(8, 2, 500, 10, i+20);
			benchmark(6, 1, 500, 10, i+30);
			benchmark(7, 1, 500, 10, i+40);
			benchmark(8, 1, 500, 5, i+50);
			benchmark(8, 1, 500, 20, i+60);		
		}
        for(int i=0; i<70; i++) {
            findMotif(i);
        }
        System.out.println("[");
        for(int i=0; i< 70; i++) {
            evaluate(i);
        }
        System.out.println("]");

    }

}