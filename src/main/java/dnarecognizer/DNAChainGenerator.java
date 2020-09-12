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
     * @param file path and name of the file that contains Original DNA chain - if you want to read from file
     * @return array of DNAChains
     */
    public ArrayList<DNAChain> Create (int DNA_Size, int OLIGONUCLEOTIDE_SIZE, int POPULATION_SIZE,StringBuilder DNA, String file){
//        try {                         // uncomment if you want to read a file
////            DNA = readFile(file);
////            System.out.println(DNA);
////
////        } catch (IOException e) {
////            e.printStackTrace();
////        }
        for(int i = 0; i<DNA_Size;i++){ // comment if you want to read a file
            switch(DNARecognizer.getGENERATOR().nextInt(4)){
                case 0:
                    DNA.append('A');
                    break;
                case 1:
                    DNA.append('C');
                    break;
                case 2:
                    DNA.append('G');
                    break;
                default:
                    DNA.append('T');
                    break;
            }
        }

        Set<String> hash_Set = new HashSet<>();
        for (int i = 0; i < DNA.length()-OLIGONUCLEOTIDE_SIZE+1; i++){
            hash_Set.add(DNA.substring(i,i+OLIGONUCLEOTIDE_SIZE));
        }
        hash_Set.remove(new StringBuilder(DNA.substring(0,OLIGONUCLEOTIDE_SIZE)));

        int k = 0;
        Iterator<String> itr = hash_Set.iterator();
        while (itr.hasNext()) {
            itr.next();
            if ((hash_Set.size()-k)>DNA_Size*(100-DNARecognizer.getPercentOfLostOligonucleotides())/100) {
                k++;
                itr.remove();
            }
            else{
                break;
            }
        }
        Set<Oligonucleotide> olig_Set = new HashSet<>();
        Oligonucleotide first = new Oligonucleotide(new StringBuilder(DNA.substring(0,OLIGONUCLEOTIDE_SIZE)));
        for(String o : hash_Set){
            Oligonucleotide oligo = new Oligonucleotide(new StringBuilder(o));
            olig_Set.add(oligo);
        }


        ArrayList<DNAChain> population = new ArrayList<>(POPULATION_SIZE);
        for(int i = 0; i < POPULATION_SIZE;i++) {
            DNAChain Chain = new DNAChain();
            ArrayList<Oligonucleotide> Olig_To_Take = new ArrayList<>(olig_Set); // list of Oligonucleotides in spectrum - all will be taken
            ArrayList<Oligonucleotide> Olig_To_Take_V2 = new ArrayList<>(olig_Set); // list of Oligonucleotides in spectrum to be used as connection to unfit ones
            ArrayList<Oligonucleotide> Olist = new ArrayList<>();
            Collections.shuffle(Olig_To_Take);
            Olist.add(0, first);
            while(Olig_To_Take.size()>0) {
                for (int j = 0; j < Olig_To_Take.size(); j++) {
                    if (Chain.fitValLoop(Olist.get(Olist.size()-1).getNucleotides(), Olig_To_Take.get(j).getNucleotides()) > 0) {
                        Olist.add(Olig_To_Take.remove(j));
                        j = -1;
                    }
                }

                if( Olig_To_Take.size() != 0){
                    for( int l = 0; l<Olig_To_Take_V2.size();l++){
                        int m1 = Chain.fitValLoop(Olist.get(Olist.size()-1).getNucleotides(),Olig_To_Take_V2.get(l).getNucleotides());
                        int m2 = Chain.fitValLoop(Olig_To_Take_V2.get(l).getNucleotides(), Olig_To_Take_V2.get(0).getNucleotides());
                        if (m1 > 0 && m2 > 2){
                            Olist.add(new Oligonucleotide(Olig_To_Take_V2.get(l).getNucleotides()));
                            Olist.add(Olig_To_Take.remove(0));
                            break;
                        }

                    }
                }


            }
            Chain.setOligonucleotides(Olist);
            population.add(Chain);
        }
        return population;
    }

}
