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
     * @param oligonucleotides the oligonucleotides to set
     */
    public void setOligonucleotides(ArrayList<Oligonucleotide> oligonucleotides) {
        this.oligonucleotides = oligonucleotides;
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
