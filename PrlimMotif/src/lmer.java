public class lmer implements Comparable<lmer>{

    int score=0;
    int c_score=0;
    StringBuilder seq;

    //seq x pos
    int[] loc;

    public lmer(String seq) {
        this.seq = new StringBuilder(seq);
    }

    public int compareTo(lmer comparelmer) {
        return comparelmer.score - this.score;
    }

}