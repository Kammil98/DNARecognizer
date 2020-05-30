package dnarecognizer.models;

public class Oligonucleotide {

    /**
     * @return the id
     */
    public int getId() {
        return id;
    }

    public Oligonucleotide(StringBuilder nucletides){
        id = idCounter++;
        this.nucleotides = nucletides;
    }
    /**
     * @return the nucleotides
     */
    public StringBuilder getNucleotides() {
        return nucleotides;
    }

    static{
        idCounter = 0;
    }

    private static int idCounter;
    private final int id;
    private final StringBuilder nucleotides;
}
