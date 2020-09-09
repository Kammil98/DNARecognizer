package dnarecognizer.models;

import java.util.ArrayList;

public class DNAChain implements Comparable<DNAChain>{

    public DNAChain(){
        id = idCounter++;
        oligonucleotides = new ArrayList<>();
    }

    /**
     * @return the id
     */
    public int getId() {
        return id;
    }

    /**
     * @return the oligonucleotides
     */
    public ArrayList<Oligonucleotide> getOligonucleotides() {
        return oligonucleotides;
    }

    /**
     *
     * @param nucl1 first nucleotide to compare
     * @param nucl2 second nucleotide to compare
     * @return the number of nucleotides that overlap between neighbour oligonucleotides
     */
    public int fitValLoop(StringBuilder nucl1, StringBuilder nucl2){
        for (int i=1; i<nucl1.length(); i++){
            for (int j=i; j<nucl1.length(); j++){
                if (nucl1.charAt(j) != nucl2.charAt(j-i)){
                    break;
                }
                else if(j == nucl1.length()-1){
                    return nucl1.length()-i;
                }

            }
        }
        return 0;
    }

    /**
     * Counts overall overlaping values between oligonucleotides in DNA chain
     */
    public void countFitVal(){
        int fitValue = 0;
        for (int i = 0; i<this.getOligonucleotides().size()-1;i++){
            fitValue+=fitValLoop(this.getOligonucleotides().get(i).getNucleotides(),this.getOligonucleotides().get(i+1).getNucleotides());
        }
        this.setFitVal(fitValue);
    }

    /**
     * @param oligonucleotides the oligonucleotides to set
     */
    public void setOligonucleotides(ArrayList<Oligonucleotide> oligonucleotides) {
        this.oligonucleotides = oligonucleotides;
        this.countFitVal();
    }

    public int getFitVal() {
        return fitVal;
    }

    public void setFitVal(int fitVal) {
        this.fitVal = fitVal;
    }

    static{
        idCounter = 0;
    }

    /**
     *
     * @return StringBuilder containing oligonucleotides joined together with overlaping
     */
    public StringBuilder toStringBuilder(){
        int f;
        StringBuilder s = new StringBuilder();
        for (int i = 0; i<this.getOligonucleotides().size()-1;i++){
            f = fitValLoop(this.getOligonucleotides().get(i).getNucleotides(),this.getOligonucleotides().get(i+1).getNucleotides());
            s.append(this.getOligonucleotides().get(i).getNucleotides().substring(0,8-f));
        }
        s.append(this.getOligonucleotides().get(this.getOligonucleotides().size()-1).getNucleotides());
        return s;
    }

    private static int idCounter;
    private final int id;
    private int fitVal = 0;
    private ArrayList<Oligonucleotide> oligonucleotides;

    @Override
    public int compareTo(DNAChain arg0) {
        return this.fitVal - arg0.getFitVal();
    }

    public void showOligonculeotidesOrder(){
        System.out.println("Oligonucleotides Order: ");
        int i = 0;
        for(Oligonucleotide olig: oligonucleotides){
            System.out.println("OligNo" + (i++) + " = "  + olig.getId());
        }
    }
}
