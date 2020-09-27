package dnarecognizer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import dnarecognizer.models.DNAChain;
import dnarecognizer.models.Oligonucleotide;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 *
 * @author kamil_2
 */
public class DNARecognizer {

    /**
     * @return the GENERATOR
     */
    public static Random getGENERATOR() {
        return GENERATOR;
    }
    public static double getPercentOfLostOligonucleotides() {return PERCENT_OF_LOST_OLIGONUCLEOTIDES;}
    public static int getOligonucleotideSize(){return OLIGONUCLEOTIDE_SIZE;}
    public static int getChildrenNoPerOneMatch(){return CHILDREN_NO_PER_ONE_MATCH;}

    //in range 300-1000 nucleotides
    private static final int DNA_SIZE = 300;
    // in range 8-10 nucleotides
    private static final int OLIGONUCLEOTIDE_SIZE =8;
    private static final int POPULATION_SIZE = 100;
    private static final int MUTATION_PROBABILITY_PERCENT = 10;
    private static final int CHILDREN_NO_PER_ONE_MATCH = 8;
    private static final float TIME_IN_MS = 30000;
    private static final double PERCENT_OF_LOST_OLIGONUCLEOTIDES = 4;
    private static final double COPIED_PERCENT_OF_DNA = 8.0d;
    private static final Random GENERATOR = new Random();
    private static final Selector SELECTOR = new Selector(POPULATION_SIZE);
    private static final DNAChainGenerator DNA_CHAIN_GENERATOR = new DNAChainGenerator();

    /**
     * selecting and crossing members of population,
     * then adding results to our population
     * @param population all possible results of our metaheuristic
     */

    public static ArrayList<DNAChain> crossover(ArrayList<DNAChain> population) {
        DNAChain member1, member2;
        ArrayList<DNAChain> group = new ArrayList<>();
        while(population.size() > 0){
            //select members to breed
            member1 = population.remove(GENERATOR.nextInt(population.size()));
            member2 = population.remove(GENERATOR.nextInt(population.size()));

            //breed
            for(int i = 0; i < CHILDREN_NO_PER_ONE_MATCH / 2; i++){
                group.add(crossTwoMembers(member1, member2, COPIED_PERCENT_OF_DNA));
                group.add(crossTwoMembers(member2, member1, COPIED_PERCENT_OF_DNA));
            }
            //add one more child, if there was odd amount of children expected
            if(CHILDREN_NO_PER_ONE_MATCH % 2 == 1)
                group.add(crossTwoMembers(member1, member2, COPIED_PERCENT_OF_DNA));

        }
        for(DNAChain D: group){
            D.countFitVal();
            population.add(D);
        }
        return population;
    }

    /**
     * crossing selected members of population.
     * @param mother main member, which will take part in crossing
     * @param father member, which will take part in crossing
     * @param copiedPartPercent how many oligs from father will be copied to childs oligs.
     *                          Given in percent
     * @return mix of Oligonucleotides order from mother and father
     */


    public static DNAChain crossTwoMembers(DNAChain mother, DNAChain father, double copiedPartPercent) {
        DNAChain child = new DNAChain();
        ArrayList<Oligonucleotide> oligs = new ArrayList<>();
        List<Oligonucleotide> motherPart1, motherPart2, motherPart3;
        int copiedPartAmount =
                (int) Math.ceil(copiedPartPercent *
                        (father.getOligonucleotides().size() / 100.0d));
        int fatherPartBegin, fatherPartEnd;

        //ensure, that ArrayList is big enough to contain all oligs (+ 5 to be certain )
        oligs.ensureCapacity(mother.getOligonucleotides().size() + 5);
        do{
            fatherPartBegin = GENERATOR.nextInt(father.getOligonucleotides().size() - copiedPartAmount - 1) + 1;}
        while(fatherPartBegin>=mother.getOligonucleotides().size());
        fatherPartEnd = fatherPartBegin + copiedPartAmount;
        if(fatherPartEnd > mother.getOligonucleotides().size() - 1){
            fatherPartEnd = mother.getOligonucleotides().size();
        }
        motherPart1 = new ArrayList<>(mother.getOligonucleotides().subList(0, fatherPartBegin));
        motherPart2 = new ArrayList<>(mother.getOligonucleotides().subList(fatherPartBegin, fatherPartEnd));
        motherPart3 = new ArrayList<>(mother.getOligonucleotides().subList(fatherPartEnd, mother.getOligonucleotides().size()));

        oligs.addAll(motherPart1);
        for(int i = 0; i <father.getOligonucleotides().size();i++){
            if(motherPart2.contains(father.getOligonucleotides().get(i))){
                oligs.add(father.getOligonucleotides().get(i));
                motherPart2.remove(father.getOligonucleotides().get(i));
            }
        }
        while(motherPart2.size()>0){
            oligs.add(motherPart2.remove(0));
        }
        oligs.addAll(motherPart3);

        child.setOligonucleotides(oligs);
        return child;
    }
    /**
     * @param population all possible results of our metaheuristic
     * @return number of best possible results in population
     */
    public static int checkLocalOptima(ArrayList<DNAChain> population){
        int optimum = population.get(0).getFitVal();
        int optimumNo = 1;
        for(DNAChain member: population){
            if(member.getFitVal() > optimum){
                optimum = member.getFitVal();
                optimumNo = 1;
            } else if(member.getFitVal() == optimum) {
                optimumNo++;
            }
        }
        return optimumNo;
    }

    /**
     * shuffle part of oligs order if local optimum is found
     * @param population all possible results of our metaheuristic
     * @param shufflingDegree how many oligs in each DNAChain should be shuffled. Given in percent
     *                        This is base value, to which we add value based on amount of similar results
     */
    public static void mutateIfLocalOptimum(ArrayList<DNAChain> population, double shufflingDegree){ //not used
        int optimumNo = checkLocalOptima(population);
        if (optimumNo < POPULATION_SIZE/20)
            return;
        Iterator<DNAChain> iter = population.iterator();
        for(DNAChain member = iter.next();
            iter.hasNext() && optimumNo > 3;
            member = iter.next()){
            mutateMemberBySwitching(member,1);
            optimumNo--;
        }
    }

    /**
     * shuffle part of oligs order in DNAChain
     * @param member DNAChain, which should be shuffled
     * @param shufflingDegree how many oligs should be shuffled. Given in percent
     */
    public static void mutateMemberByShuffling(DNAChain member, double shufflingDegree) {
        ArrayList<Oligonucleotide> oliginucleotides = member.getOligonucleotides();
        int shufflingAmount =
                (int) Math.ceil(shufflingDegree *
                        (oliginucleotides.size() / 100.0d));
        int shufflingBegin = GENERATOR.nextInt(
                oliginucleotides.size() - shufflingAmount -1 ) + 1;
        List<Oligonucleotide> oligToShuffle = oliginucleotides
                .subList(shufflingBegin, shufflingBegin + shufflingAmount);
        Collections.shuffle(oligToShuffle);
    }

    /**
     * switches places of random oligs in DNAChain
     * @param member DNAChain, which should have its oligs switched
     * @param shufflingAmount how many oligs to switch places
     */
    public static void mutateMemberBySwitching(DNAChain member, int shufflingAmount) {
        ArrayList<Oligonucleotide> oliginucleotides = member.getOligonucleotides();
        for(int i = 0;i<shufflingAmount;i++) {
            int[] switcho = new int[2];
            for(int s = 0; s<switcho.length;s++) {
                do {
                    do {
                        switcho[s] = GENERATOR.nextInt(oliginucleotides.size() - 1) + 1;
                    } while (member.fitValLoop(oliginucleotides.get(switcho[s] - 1).getNucleotides(), oliginucleotides.get(switcho[s]).getNucleotides()) == OLIGONUCLEOTIDE_SIZE - 1);
                    if (switcho[s] + 1 == member.getOligonucleotides().size()) {
                        break;
                    }
                } while (member.fitValLoop(oliginucleotides.get(switcho[s]).getNucleotides(), oliginucleotides.get(switcho[s] + 1).getNucleotides()) == OLIGONUCLEOTIDE_SIZE - 1||switcho[0]==switcho[1]);
            }
            Collections.swap(oliginucleotides,switcho[0],switcho[1]);
        }
    }

    public static void mutateMemberByChangingPosition(DNAChain member) {
        ArrayList<Oligonucleotide> oliginucleotides = member.getOligonucleotides();
        int changed, place;
        do {
            do {
                changed = GENERATOR.nextInt(oliginucleotides.size() - 1) + 1;
            } while (member.fitValLoop(oliginucleotides.get(changed - 1).getNucleotides(), oliginucleotides.get(changed).getNucleotides()) == OLIGONUCLEOTIDE_SIZE - 1);
            if (changed + 1 == member.getOligonucleotides().size()) {
                break;
            }
        } while (member.fitValLoop(oliginucleotides.get(changed).getNucleotides(), oliginucleotides.get(changed + 1).getNucleotides()) == OLIGONUCLEOTIDE_SIZE - 1);
        do{
            place = GENERATOR.nextInt(oliginucleotides.size() - 1) + 1;
        }  while(member.fitValLoop(oliginucleotides.get(place - 1).getNucleotides(), oliginucleotides.get(place).getNucleotides())== OLIGONUCLEOTIDE_SIZE - 1);
        if(changed>place){
            oliginucleotides.add(place,oliginucleotides.remove(changed));
        }else{
            oliginucleotides.add(place - 1,oliginucleotides.remove(changed));
        }
    }

    public static void mutateMemberByChangingSectionPosition(DNAChain member) {
        ArrayList<Oligonucleotide> oliginucleotides = member.getOligonucleotides();
        for(int i = 1;i<oliginucleotides.size()-1;i++){
            if(member.fitValLoop(oliginucleotides.get(i).getNucleotides(), oliginucleotides.get(i+1).getNucleotides()) == OLIGONUCLEOTIDE_SIZE-1){
                break;
            }
            else if (i+1==oliginucleotides.size()-1) {
                return;
            }
        }
        int change1, change2, position;
        do {
            change1 = GENERATOR.nextInt(oliginucleotides.size() - 2) + 1;
        }while (member.fitValLoop(oliginucleotides.get(change1 - 1).getNucleotides(), oliginucleotides.get(change1).getNucleotides()) != OLIGONUCLEOTIDE_SIZE-1 &&
                member.fitValLoop(oliginucleotides.get(change1).getNucleotides(), oliginucleotides.get(change1 + 1).getNucleotides()) != OLIGONUCLEOTIDE_SIZE-1);

        change2 = change1;

        while(member.fitValLoop(oliginucleotides.get(change1-1).getNucleotides(),oliginucleotides.get(change1).getNucleotides())==OLIGONUCLEOTIDE_SIZE-1){
            change1--;
            if(change1<1){
                return;
            }
        }
        if(change1 == 0){
            System.out.println("------------------------------------------------------------------------------------------");
        }

        while (member.fitValLoop(oliginucleotides.get(change2).getNucleotides(),oliginucleotides.get(change2+1).getNucleotides())==OLIGONUCLEOTIDE_SIZE-1){
            change2++;
            if(change2>=oliginucleotides.size()-1){
                break;
            }
        }

        do{
            position = GENERATOR.nextInt(oliginucleotides.size() - 1) + 1;
        }while(!(position < change1 || position > change2) ||member.fitValLoop(oliginucleotides.get(position-1).getNucleotides(),oliginucleotides.get(position).getNucleotides())==OLIGONUCLEOTIDE_SIZE-1);
        if(position == 0){
            System.out.println("p------------------------------------------------------------------------------------------");
        }
        ArrayList<Oligonucleotide> newchain = new ArrayList<>();

        if(position<change1) {
            newchain.addAll(oliginucleotides.subList(0, position));
            newchain.addAll(oliginucleotides.subList(change1, change2+1));
            newchain.addAll(oliginucleotides.subList(position, change1));
            newchain.addAll(oliginucleotides.subList(change2+1, oliginucleotides.size()));
        }else{
            newchain.addAll(oliginucleotides.subList(0, change1));
            newchain.addAll(oliginucleotides.subList(change2+1,position));
            newchain.addAll(oliginucleotides.subList(change1,change2+1));
            newchain.addAll(oliginucleotides.subList(position, oliginucleotides.size()));
        }
        member.setOligonucleotides(newchain);

    }


    /**
     * Add some oligs to result
     * @param member result, which will by mutated
     * @param mutationSize how many oligonucleotides should be mutated. Given in percent
     */
    public static void mutateMemberByAdding(DNAChain member, double mutationSize){
        int mutationLocation, indexOfLastChar;
        int mutationAmount =
                (int) Math.ceil(mutationSize *
                        (member.getOligonucleotides().size() / 100.0d));
        StringBuilder oligVal;
        char first, last, nucleotide;
        for(int mutationNo = 0; mutationNo < mutationAmount; mutationNo++){
            mutationLocation = 1 + GENERATOR.nextInt(OLIGONUCLEOTIDE_SIZE - 1);
            //make sure first and last char in oligonculeotide
            //will be the same, as in oligonculeotides next to it
            indexOfLastChar = member.getOligonucleotides()
                    .get(mutationLocation - 1)
                    .getNucleotides().length() - 1;
            first = member.getOligonucleotides()
                    .get(mutationLocation - 1)
                    .getNucleotides().charAt(indexOfLastChar);
            last = member.getOligonucleotides()
                    .get(mutationLocation)
                    .getNucleotides().charAt(0);

            //Create oligonculeotide
            oligVal = new StringBuilder().append(first);
            for(int oligNo = 1; oligNo < OLIGONUCLEOTIDE_SIZE - 1; oligNo++){
                switch(GENERATOR.nextInt(4)){
                    case 0:
                        nucleotide = 'A';
                        break;
                    case 1:
                        nucleotide = 'C';
                        break;
                    case 2:
                        nucleotide = 'G';
                        break;
                    default:
                        nucleotide = 'T';
                        break;
                }
                oligVal.append(nucleotide);
            }
            oligVal.append(last);
            member.getOligonucleotides().add(mutationLocation,
                    new Oligonucleotide(oligVal));
        }
    }

    /**
     *
     * @param population all possible results of our metaheuristic
     * @param progress elapsed time given in percent
     */
    public static void select(ArrayList<DNAChain> population, float progress){
        if(progress < 99){
            SELECTOR.contest(population);
        }
        else{
            SELECTOR.ranking(population);
        }
    }

    /**
     * for now prints fit value of first DNA chain in curent population and at the end shows best (for fitValue) DNA chain its lenght, all nucleodies in string and all oligonucleotides separatly.
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        int AvgBuildLenght = 0, AvgLevis = 0, AvgOligNum = 0, AvgFitVal = 0;
        for (int z = 0; z < 10; z++) {

            long startTime = System.currentTimeMillis();
            int optnr = 0, optval = 0; // optnr - number of occurrences of the same best fitval, optval - fitval of the best chain in previous loop occurrence
            DNAChain bestchain = new DNAChain(); // best chain in generation
            DNAChain supreme = new DNAChain(); // best chain in all generations
            ArrayList<DNAChain> population;
            StringBuilder DNA = new StringBuilder();
            population = DNA_CHAIN_GENERATOR.Create(DNA_SIZE, OLIGONUCLEOTIDE_SIZE, POPULATION_SIZE, DNA, "src/main/resources/SampleDNA.txt");
            while ((System.currentTimeMillis() - startTime) < TIME_IN_MS) {

                for (DNAChain dnaChain : population) {
                    if (dnaChain.getFitVal() > bestchain.getFitVal()) {
                        bestchain = dnaChain;
                    }
                }

                if (bestchain.getFitVal() > supreme.getFitVal()) {
//                    System.out.println("New supreme: " + bestchain.getFitVal());
//                    System.out.println("Time taken: " + (System.currentTimeMillis() - startTime) + "ms");
                    supreme.setOligonucleotides(bestchain.getOligonucleotides());
                    optnr = 0;
                }
                if (bestchain.getFitVal() == optval)
                    optnr++;
                else if (bestchain.getFitVal() > optval)
                    optval = bestchain.getFitVal();

                if (optnr >= 30) {
                    for (DNAChain dnaChain : population) {
                        if (dnaChain.getOligonucleotides().size() < DNA_SIZE * 1.01)
                            mutateMemberByAdding(dnaChain, 1);
                        else

                            for (int l = 0; l < 5; l++) {
                                mutateMemberByChangingPosition(dnaChain);
                                if (((System.currentTimeMillis() - startTime) / TIME_IN_MS * 100) > 50)
                                    mutateMemberByChangingSectionPosition(dnaChain);
                            }
                    }

                    for (DNAChain D : population) {
                        D.countFitVal();
                    }
                    optnr = 0;
                    optval = 0;
                    bestchain = new DNAChain();
                }
                crossover(population);
                for (DNAChain dnaChain : population) {
                    if (GENERATOR.nextInt(100) < MUTATION_PROBABILITY_PERCENT ) {
                        mutateMemberByChangingPosition(dnaChain);
                        if (((System.currentTimeMillis() - startTime) / TIME_IN_MS * 100) > 50) {
                            mutateMemberByChangingSectionPosition(dnaChain);
                        }
                    }
                }
                for (DNAChain D : population) {
                    D.countFitVal();
                }
                select(population, (System.currentTimeMillis() - startTime) / TIME_IN_MS * 100);
                population.remove(GENERATOR.nextInt(population.size()));
                DNAChain OG = new DNAChain();
                OG.setOligonucleotides(supreme.getOligonucleotides());
                OG.countFitVal();
                population.add(OG);


            }
            for (int i = 0; i < supreme.getOligonucleotides().size(); i++) {
                if (i != supreme.getOligonucleotides().size() - 1) {
                    System.out.println("Olig #" + i + ": " + supreme.getOligonucleotides().get(i).getNucleotides() + " --> f: " + supreme.fitValLoop(supreme.getOligonucleotides().get(i).getNucleotides(), supreme.getOligonucleotides().get(i + 1).getNucleotides()));
                } else {
                    System.out.println("Olig #" + i + ": " + supreme.getOligonucleotides().get(i).getNucleotides());
                }
            }
            System.out.println("Best fitValue " + supreme.getFitVal());
            System.out.println("How many Oligs: " + supreme.getOligonucleotides().size());
            System.out.println("Build DNA lenght: " + supreme.toStringBuilder().length());
            System.out.println("Build DNA: " + supreme.toStringBuilder());
            System.out.println("Orig. DNA: " + DNA);
            System.out.println("Time taken: " + (System.currentTimeMillis() - startTime) + "ms");
            System.out.println("LevenshteinDistance: " + LevenshteinDistance.Compute(supreme.toString(), DNA.toString()));
            AvgFitVal+=supreme.getFitVal();
            AvgOligNum+=supreme.getOligonucleotides().size();
            AvgBuildLenght+=supreme.toStringBuilder().length();
            AvgLevis+=LevenshteinDistance.Compute(supreme.toString(), DNA.toString());

            try {
                File myObj = new File("Test_Mutation_pos_"+MUTATION_PROBABILITY_PERCENT+"_"+ z+".txt");
                if (myObj.createNewFile()) {
                    System.out.println("File created: " + myObj.getName());
                } else {
                    System.out.println("File already exists.");
                }
            } catch (IOException e) {
                System.out.println("An error occurred in creating");
                e.printStackTrace();
            }
            try {
                FileWriter myWriter = new FileWriter("Test_Mutation_pos_"+MUTATION_PROBABILITY_PERCENT+"_"+ z+".txt");

                myWriter.write("Best fitValue " + supreme.getFitVal()+"\n");
                myWriter.write("How many Oligs: " + supreme.getOligonucleotides().size()+"\n");
                myWriter.write("Build DNA lenght: " + supreme.toStringBuilder().length()+"\n");
                myWriter.write("Build DNA: " + supreme.toStringBuilder()+"\n");
                myWriter.write("Orig. DNA: " + DNA+"\n");
                myWriter.write("Time taken: " + (System.currentTimeMillis() - startTime) + "ms"+"\n");
                myWriter.write("LevenshteinDistance: " + LevenshteinDistance.Compute(supreme.toString(), DNA.toString())+"\n");
                myWriter.close();
                System.out.println("Successfully wrote to the file.");
            } catch (IOException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
        }
        try {
            File myObj = new File("Test_Mutation_pos_"+MUTATION_PROBABILITY_PERCENT+"_Total.txt");
            if (myObj.createNewFile()) {
                System.out.println("File created: " + myObj.getName());
            } else {
                System.out.println("File already exists.");
            }
        } catch (IOException e) {
            System.out.println("An error occurred in creating");
            e.printStackTrace();
        }
        try {
            FileWriter myWriter = new FileWriter("Test_Mutation_pos_"+MUTATION_PROBABILITY_PERCENT+"_Total.txt");
            myWriter.write("AvgFitVal: "+AvgFitVal/10+"\n");
            myWriter.write("AvgOligNum: "+AvgOligNum/10+"\n");
            myWriter.write("AvgBuildLenght: "+AvgBuildLenght/10+"\n");
            myWriter.write("AvgLevis: "+AvgLevis/10+"\n");
            myWriter.close();
            System.out.println("Successfully wrote to the file.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

    }

}
