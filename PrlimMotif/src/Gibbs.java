import java.io.*;
import java.util.LinkedList;
import java.util.List;

public class Gibbs {

    public static void main(String[] args) {
        //ml, nm, sl, sc
		/*
		for(int i=0; i<10; i++)
		{
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
            findMotif( 16 );
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void findMotif(int setnum) throws IOException {
        // Will try to implement with Gibbs Sampling as described here: http://bix.ucsd.edu/bioalgorithms/presentations/Ch12_RandAlgs.pdf

        String[] sequencesStrings = getSequences(setnum);
        int motifLength           = getMotifLength(setnum);

        Sequence[] sequences = new Sequence[sequencesStrings.length];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = new Sequence(motifLength, sequencesStrings[i]);
        }

        float[][] oldprofile = new float[motifLength][4];
        float[][] profile    = new float[motifLength][4];
        fillArray(profile,0f);

        int differenceThresh = 1, newscore = 0, oldscore = Integer.MIN_VALUE + differenceThresh;
        int iterations = 0;
        while ( iterations < 90000) {

            iterations++;
            //System.out.println("---------------------------------\ni : " + iterations + "\nNewScore : " + newscore + "\nOldScore : " + oldscore);
            int removedSequenceIndex = (int) Math.floor( Math.random() * sequences.length );

            // Fill Profile Matrix
            for (int sequenceIndex = 0; sequenceIndex < sequences.length; sequenceIndex++) {

                if (sequenceIndex == removedSequenceIndex)
                    continue;

                Sequence seq     = sequences[sequenceIndex];
                String currMotif = seq.getCurrentMotif();

                for (int motifIndex = 0; motifIndex < motifLength; motifIndex++) {

                    char base = currMotif.charAt( motifIndex );
                    profile[ motifIndex ][ baseToInt( base ) ] += ( 1f / (float) sequences.length );

                }

            }

            //printArray(profile);

            // Calculate the prob(a|P) for every possible 8-mer in the removed sequence ... a = sequences[removedSequence]
            float[] probDist = new float[ sequences[ removedSequenceIndex ].possibleMotifs.length ];
            for (int i = 0; i < probDist.length; i++)
                probDist[i] = 0f; // Initialize all values to 0

            for (int i = 0; i < sequences[removedSequenceIndex].possibleMotifs.length; i++) {

                String motif = sequences[removedSequenceIndex].possibleMotifs[i];

                // create probability distribution for next starting point in sequences[removedSequenceIndex]
                for (int character = 0; character < motif.length(); character++) {
                    /*
                    System.out.println("Profile \'A\' : " + profile[character][0]);
                    System.out.println("Profile \'C\' : " + profile[character][1]);
                    System.out.println("Profile \'G\' : " + profile[character][2]);
                    System.out.println("Profile \'T\' : " + profile[character][3]);
                    */
                    switch ( motif.charAt(character) ) {
                        case 'A':
                            probDist[i] += profile[character][0] == 0 ? 0 : Math.log(profile[character][0] / (sequences[removedSequenceIndex].Afreq / sequences[removedSequenceIndex].length() ) ); // 0 == A
                            break;
                        case 'C':
                            probDist[i] += profile[character][1] == 0 ? 0 : Math.log(profile[character][1] / (sequences[removedSequenceIndex].Cfreq / sequences[removedSequenceIndex].length() ) ); // 1 == C
                            break;
                        case 'G':
                            probDist[i] += profile[character][2] == 0 ? 0 : Math.log(profile[character][2] / (sequences[removedSequenceIndex].Gfreq / sequences[removedSequenceIndex].length() ) ); // 2 == G
                            break;
                        case 'T':
                            probDist[i] += profile[character][3] == 0 ? 0 : Math.log(profile[character][3] / (sequences[removedSequenceIndex].Tfreq / sequences[removedSequenceIndex].length() ) ); // 3 == T
                            break;
                        default:
                            System.out.println("We got anotha problem");
                            break;
                    }

                }

            }

            // Find lowest value in probDist
            float lowest = probDist[0];

            for ( float value : probDist ) {
                if (value < lowest)
                    lowest = value;
            }

            //System.out.println("Lowest : " + lowest);

            // Sort of normalizing the probability distribution so everything is greater than or equal to 1
            float sum = 0;
            for (int i = 0; i < probDist.length; i++) {
                probDist[i] = probDist[i]/lowest;
                sum += probDist[i];
            }

            for (int i = 0; i < probDist.length; i++) {
                probDist[i] /= sum;
            }


            float probValue = (float) ( Math.random() * 1f ); // chooses a new starting position for the sequence based on the probability distribution
            float startProb = 0, endProb = probDist[0];
            for (int i = 0; i < probDist.length - 1; i++) {

                if ( probValue >= startProb && probValue <= endProb ) {
                    // reset the motifStarting point to i
                    sequences[removedSequenceIndex].setMotifStartIndex(i);
                    break;
                }

                startProb = startProb + probDist[i];
                endProb   = endProb + probDist[i+1];
                //System.out.println("startProb : " + startProb + " , endProb : " + endProb);
            }

            // Calculate new and oldscores
            oldscore = newscore;
            newscore = 0;
            for (int x = 0; x < profile.length; x++) {
                for (int y = 0; y < profile[0].length; y++) {
                    newscore += profile[x][y] - oldprofile[x][y];
                }
            }
            copyArray(profile,oldprofile);
            fillArray(profile,0f);
        }

        System.out.println("Iterations : " + iterations);

        for (Sequence seq : sequences) {
            System.out.println(seq.getCurrentMotif().replaceAll("[.]*","  "));
        }

    }

    public static int baseToInt(char base) {
        switch (base) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                System.out.println("Oh No!");
                return -1;
        }
    }

    public static char intToBase(int i) {
        switch (i) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                System.out.println(":'(");
                return '\0';
        }
    }

    public static void copyArray(float[][] from, float[][] to) {
        for (int x = 0; x < from.length; x++) {
            for (int y = 0; y < from[0].length; y++) {
                to[x][y] = from[x][y];
            }
        }
    }

    public static void fillArray(float[][] grid, float fill) {
        for (int x = 0; x < grid.length; x++) {
            for (int y = 0; y < grid[0].length; y++) {
                grid[x][y] = fill;
            }
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
        private final String[] possibleMotifs;
        public final int Afreq;
        public final int Cfreq;
        public final int Gfreq;
        public final int Tfreq;

        public Sequence(int motifLength, String sequence) {
            this.sequence = sequence;
            this.possibleMotifs = new String[sequence.length() - motifLength];
            // fill in this.possibleMotifs with all possible l-mers of length == motifLength
            for (int start =  0, end = motifLength; end < sequence.length();) {

                this.possibleMotifs[start] =  this.sequence.substring(start, end);

                start++;
                end++;
            }

            this.motifStartIndex = (int) Math.floor( Math.random() * this.possibleMotifs.length );

            // find out amino acid frequencies
            int[] freqs = new int[4];
            for (int i = 0; i < this.sequence.length(); i++) {
                freqs[ baseToInt( this.sequence.charAt(i) ) ]++;
            }
            this.Afreq = freqs[0];
            this.Cfreq = freqs[1];
            this.Gfreq = freqs[2];
            this.Tfreq = freqs[3];

        }

        public int length() {
            return this.sequence.length();
        }

        public String getCurrentMotif() {
            return this.possibleMotifs[this.motifStartIndex];
        }

        public char charAt(int index) {
            return this.sequence.charAt(index);
        }

        public void setMotifStartIndex(int newStart) {
            this.motifStartIndex = newStart;
        }


    }

    /**
     * Pretty prints a 2D Number Array to console ... For debugging only
     * @param array
     */
    private static void printArray( float[][] array ) {
        int maxNumStringLength = 9;
        for (int x = 0; x < array.length; x++) {
            String str = "";
            for (int y = 0; y < array[0].length; y++) {
                String num = String.valueOf(array[x][y]);
                for (int i = 0; i < maxNumStringLength - num.length(); i++) {
                    num = " " + num;
                }
                str += num + " ";
            }
            System.out.println(str);
        }
    }

}
