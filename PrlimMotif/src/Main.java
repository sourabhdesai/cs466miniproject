import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


public class Main {

    public static void main(String[] args) {
        //ml, nm, sl, sc
		/*
		for(int i=0; i<10; i++) {
			benchmark(8, 1, 500, 10, i);
			benchmark(8, 0, 500, 10, i+10);
			benchmark(8, 2, 500, 10, i+20);
			benchmark(6, 1, 500, 10, i+30);
			benchmark(7, 1, 500, 10, i+40);
			benchmark(8, 1, 500, 5, i+50);
			benchmark(8, 1, 500, 20, i+60);
		}
	    */

        try {
            findMotif( 0 );
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

	public static void findMotif(int setnum) throws IOException {
        // Will try to implement with Gibbs Sampling as described here: http://bix.ucsd.edu/bioalgorithms/presentations/Ch12_RandAlgs.pdf

        String[] sequencesStrings = getSequences(setnum);
        int motifLength = getMotifLength(setnum);

        Sequence[] sequences = new Sequence[sequencesStrings.length];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = new Sequence(motifLength, sequencesStrings[i]);
        }

        float[][] profile = new float[motifLength][4];

        int[] startPoints = new int[sequences.length];

        while (true) {

            for (int i = 0; i < sequences.length; i++) {
                startPoints[i] = (int) Math.floor( Math.random() * sequences[i].length() );
            }

            int removedSequence = (int) Math.floor( Math.random() * sequences.length );

            // Fill Profile Matrix
            for (int sequenceIndex = 0; sequenceIndex < sequences.length; sequenceIndex++) {

                if (sequenceIndex == removedSequence)
                    continue;

                Sequence seq     = sequences[sequenceIndex];
                int startPoint = startPoints[sequenceIndex];

                for (int motifIndex = 0; motifIndex < motifLength; motifIndex++) {

                    char base = seq.charAt( motifIndex + startPoint );

                    switch (base) {
                        case 'A':
                            profile[motifIndex][0] += (float) ( 1 / sequences.length );
                            break;
                        case 'C':
                            profile[motifIndex][1] += (float) ( 1 / sequences.length );;
                            break;
                        case 'G':
                            profile[motifIndex][2] += (float) ( 1 / sequences.length );
                            break;
                        case 'T':
                            profile[motifIndex][3] += (float) ( 1 / sequences.length );
                            break;
                        default:
                            System.out.println("Houston we have a problem...");
                            break;
                    }

                }

            }

            // TODO: Rest of Gibbs Sampling Algorithm ... Refer to link above for help


        }
	}

    public static String[] getSequences(int setnum) throws IOException {
        File seqF = new File( "data/set" + setnum + "/sequences.fa" );

        BufferedReader seqR =  new BufferedReader( new FileReader(seqF) );

        List<String> sequences = new LinkedList<String>();
        String currentLine;

        while( ( currentLine = seqR.readLine() ) != null ) {
            if (currentLine.charAt(0) != '>')
                sequences.add(currentLine);
        }

        seqR.close();

        return sequences.toArray( new String[ sequences.size() ] );
    }

    public static int getMotifLength(int setnum) throws IOException {
        File mlenF = new File("data/set"+setnum+"/motiflength.txt");
        BufferedReader mlenR =  new BufferedReader(new FileReader(mlenF));
        int mlen = Integer.parseInt(mlenR.readLine());
        mlenR.close();
        return mlen;
    }

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

    public static class Sequence {

        private int motifStartIndex;
        private final String sequence;
        private final List<String> possibleMotifs;

        public Sequence(int motifLength, String sequence) {
            this.sequence = sequence;
            this.possibleMotifs = new ArrayList<String>();
            // fill in this.possibleMotifs with all possible l-mers of length == motifLength
            for (int start =  0, end = motifLength; end < sequence.length();) {

                this.possibleMotifs.add( this.sequence.substring(start, end) );

                start++;
                end++;
            }

            this.motifStartIndex = (int) Math.floor( Math.random() * this.sequence.length() );

        }

        public int length() {
            return this.sequence.length();
        }

        public String getCurrentMotif() {
            return this.possibleMotifs.get(this.motifStartIndex);
        }

        public char charAt(int index) {
            return this.sequence.charAt(index);
        }


    }

}
