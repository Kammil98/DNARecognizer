package dnarecognizer;

import dnarecognizer.models.DNAChain;
import dnarecognizer.models.Oligonucleotide;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;


public class DNAChainGenerator {
    /**
     *
     * @param file path and name to the file that contains Original DNA chain
     * @return StringBuilder containing read nucleotides
     */
    public StringBuilder readFile(String file) throws IOException {

        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            StringBuilder stringBuilder = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                stringBuilder.append(line);
            }

            return stringBuilder;
        }
    }

    /**
     *
     * @param DNA_Size number of nucleotides in DNA chain
     * @param OLIGONUCLEOTIDE_SIZE number of nucleotides in oligonucleotide
     * @param POPULATION_SIZE number of DNA chain instances
     * @param file path and name of the file that contains Original DNA chain
     * @return array of DNAChains
     */
    public ArrayList<DNAChain> Create (int DNA_Size, int OLIGONUCLEOTIDE_SIZE, int POPULATION_SIZE, String file){
        StringBuilder DNA = new StringBuilder();
        try {
            DNA = readFile(file);
            System.out.println(DNA);

        } catch (IOException e) {
            e.printStackTrace();
        }
        Set<String> hash_Set = new HashSet<>();
        for (int i = 0; i < DNA.length()-OLIGONUCLEOTIDE_SIZE+1; i++){
            hash_Set.add(DNA.substring(i,i+OLIGONUCLEOTIDE_SIZE));
        }
        System.out.println("-----------------  Hashsize: "+hash_Set.size());
        hash_Set.remove(new StringBuilder(DNA.substring(0,OLIGONUCLEOTIDE_SIZE)));

        Set<Oligonucleotide> olig_Set = new HashSet<>();
        Oligonucleotide first = new Oligonucleotide(new StringBuilder(DNA.substring(0,OLIGONUCLEOTIDE_SIZE)));
        for(String o : hash_Set){
            Oligonucleotide oligo = new Oligonucleotide(new StringBuilder(o));
            olig_Set.add(oligo);
        }
        ArrayList<DNAChain> population = new ArrayList<>(POPULATION_SIZE);
        for(int i = 0; i < POPULATION_SIZE;i++){
            DNAChain Chain = new DNAChain();
            ArrayList<Oligonucleotide> Olist = new ArrayList<>(olig_Set);
            Collections.shuffle(Olist);
            Olist.add(0,first);
            Chain.setOligonucleotides(Olist);
//            for(int x = 0; x  < 5;x++)
//                System.out.println(Chain.getOligonucleotides().get(x).getNucleotides());
            population.add(Chain);
        }
        return population;
    }

}