package dnarecognizer;

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

    //in range 300-1000 nucleotides
    private static final int DNA_SIZE = 300;
    // in range 8-10 nucleotides
    private static final int OLIGONUCLEOTIDE_SIZE = 8;
    private static final int POPULATION_SIZE = 100;
    private static final int EXPANDED_POPULATION_SIZE = 400;
    private static final int NUMBER_OF_GENERATIONS = 2000;
    private static final int CHILDREN_NO_PER_ONE_MATCH = 4;
    private static final double PERCENT_OF_LOST_OLIGONUCLEOTIDES = 4;
    private static final double COPIED_PERCENT_OF_DNA = 10.0d;
    private static final Random GENERATOR = new Random();
    private static final Selector SELECTOR = new Selector(POPULATION_SIZE);
    private static final DNAChainGenerator DNA_CHAIN_GENERATOR = new DNAChainGenerator();

    /**
     * selecting and crossing members of population,
     * then adding results to our population
     * @param population all possible results of our metaheuristic
     */
    public static ArrayList<DNAChain> crossover(ArrayList<DNAChain> population) {
        int breedNo = EXPANDED_POPULATION_SIZE - population.size();
        int groupSize = POPULATION_SIZE / 20;
        DNAChain member1, member2;
        while(breedNo > 0){
            ArrayList<DNAChain> group = new ArrayList<>();

            //select members to breed
            for(int i = 0; i < groupSize; i++)
                group.add(population.get(
                        GENERATOR.nextInt(population.size())));
            member1 = SELECTOR.bestOfGroup(group);
            group.remove(member1);
            member2 = SELECTOR.bestOfGroup(group);

            //breed
            for(int i = 0; i < CHILDREN_NO_PER_ONE_MATCH / 2; i++){
                population.add(crossTwoMembers(member1, member2, COPIED_PERCENT_OF_DNA));
                population.add(crossTwoMembers(member2, member1, COPIED_PERCENT_OF_DNA));
            }
            //add one more child, if there was odd amount of children expected
            if(CHILDREN_NO_PER_ONE_MATCH % 2 == 1)
                population.add(crossTwoMembers(member1, member2, COPIED_PERCENT_OF_DNA));

            breedNo -= CHILDREN_NO_PER_ONE_MATCH;
        }
        for(DNAChain D: population){
            D.countFitVal();
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
        List<Oligonucleotide> fatherPart2, motherPart1, motherPart2, motherPart3;
        int copiedPartAmount =
                (int) Math.ceil(copiedPartPercent *
                        (father.getOligonucleotides().size() / 100.0d));
        int fatherPartBegin, fatherPartEnd;

        //ensure, that ArrayList is big enough to contain all oligs (+ 5 to be certain )
        oligs.ensureCapacity(mother.getOligonucleotides().size() + 5);

        fatherPartBegin = GENERATOR.nextInt(father.getOligonucleotides().size() - copiedPartAmount - 1) + 1;
        fatherPartEnd = fatherPartBegin + copiedPartAmount;
        if(fatherPartEnd > mother.getOligonucleotides().size() - 1){
            fatherPartEnd = mother.getOligonucleotides().size();
        }
        motherPart1 = new ArrayList<>(mother.getOligonucleotides().subList(0, fatherPartBegin));
        motherPart2 = new ArrayList<>(mother.getOligonucleotides().subList(fatherPartBegin, fatherPartEnd));
        motherPart3 = new ArrayList<>(mother.getOligonucleotides().subList(fatherPartEnd, mother.getOligonucleotides().size()));
        fatherPart2 = new ArrayList<>(father.getOligonucleotides().subList(fatherPartBegin, fatherPartEnd));

        //remove from father part all oligs, which already appeared in mothers part 1 or 3
        fatherPart2.removeAll(motherPart1);
        fatherPart2.removeAll(motherPart3);

        //left only this oligs, which didn't appeared in fathers second part
        motherPart2.removeAll(fatherPart2);

        //shallow copy - the same oligs in Parents and child (but another ArrayLists of oligs
        oligs.addAll(motherPart1);
        oligs.addAll(motherPart2);
        oligs.addAll(fatherPart2);
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
    public static void mutateIfLocalOptimum(ArrayList<DNAChain> population,
                                            double shufflingDegree){
        int optimumNo = checkLocalOptima(population);
        //don't need to
        if (optimumNo < POPULATION_SIZE/20)
            return;
        Iterator<DNAChain> iter = population.iterator();
        for(DNAChain member = iter.next();
            iter.hasNext() && optimumNo > 3;
            member = iter.next()){
//            if(optimumNo > POPULATION_SIZE / 5)
//                shufflingDegree = Math.min(5.0d, shufflingDegree + 2.0d); // changed from max to min cause if 100.0d then its always 100% of oligonucleotides mutated
//            else                                                            // changed numbers from 100.0d and 70.0d to 5.0d and 2.0d
//                shufflingDegree = Math.min(3.0d, shufflingDegree + 0.0d);  // changed numbers from 100.0d and 30.0d to 3.0d and 0.0d
            //mutateMemberByShuffling(member, shufflingDegree);
            mutateMemberBySwitching(member,1);
            optimumNo--;
        }
//        System.out.println("Shuffle Kamil's");
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
    public static void mutateMemberBySwitching(DNAChain member, double shufflingAmount) {
        ArrayList<Oligonucleotide> oliginucleotides = member.getOligonucleotides();
        for(int i = 0;i<shufflingAmount;i++){
            int switch1 = GENERATOR.nextInt(oliginucleotides.size() -1 ) + 1;
            int switch2;
            do {
                switch2 = GENERATOR.nextInt(oliginucleotides.size() - 1) + 1;
            }while  (member.fitValLoop(oliginucleotides.get(switch1).getNucleotides(), oliginucleotides.get(switch1 - 1).getNucleotides())==0&&
                    member.fitValLoop(oliginucleotides.get(switch2).getNucleotides(), oliginucleotides.get(switch2 - 1).getNucleotides())==0
                    );
            Collections.swap(oliginucleotides,switch1,switch2);
        }
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
        if(progress < 5){
            SELECTOR.roulette(population);
        }
        else if(progress <= 100){
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
        long startTime = System.currentTimeMillis();
        int optnr = 0, optval = 0; // optnr - number of occurrences of the same best fitval, optval - fitval of the best chain in previous loop occurrence
        DNAChain bestchain = new DNAChain(); // best chain in generation
        DNAChain supreme = new DNAChain(); // best chain in all generations
        ArrayList<DNAChain> population;
        StringBuilder DNA = new StringBuilder();
        population = DNA_CHAIN_GENERATOR.Create(DNA_SIZE,OLIGONUCLEOTIDE_SIZE,POPULATION_SIZE,DNA,"src/main/resources/SampleDNA.txt");

        for (int i=0;i<NUMBER_OF_GENERATIONS;i++){
            System.out.println("NR: "+i);
            for (DNAChain dnaChain : population) {
                if (dnaChain.getFitVal() > bestchain.getFitVal()) {
                    bestchain = dnaChain;
                }
            }
            System.out.println("Fit Value gen. best: " + bestchain.getFitVal());
            if(bestchain.getFitVal()>supreme.getFitVal()){
                System.out.println("New supreme: " + bestchain.getFitVal());
                supreme.setOligonucleotides(bestchain.getOligonucleotides());
                optnr=0;
            }
            if(bestchain.getFitVal() == optval)
                optnr++;
            else if(bestchain.getFitVal() > optval)
                optval = bestchain.getFitVal();

            if(optnr>=30){
                for (DNAChain dnaChain : population) {
                    if(population.size()<DNA_SIZE*0.01)
                    mutateMemberByAdding(dnaChain,1);
                    else
                    mutateMemberBySwitching(dnaChain, 5);
                }
//                System.out.println("Shuffle Filip's");
                for(DNAChain D: population){
                    D.countFitVal();
                }
                optnr = 0;
                optval = 0;
                bestchain = new DNAChain();
            }
            crossover(population);
            select(population,i/NUMBER_OF_GENERATIONS*100);
            mutateIfLocalOptimum(population,3);

            for(DNAChain D: population){
                D.countFitVal();
            }

        }
        for(int i = 0; i<supreme.getOligonucleotides().size();i++){
            System.out.println(supreme.getOligonucleotides().get(i).getNucleotides());
        }
        System.out.println("Best fitValue " + supreme.getFitVal());
        System.out.println("How many Oligs: " + supreme.getOligonucleotides().size());
        System.out.println("Build DNA lenght: " + supreme.toStringBuilder().length());
        System.out.println("Build DNA: " + supreme.toStringBuilder());
        System.out.println("Orig. DNA: " + DNA);
        System.out.println("Time taken: "+(System.currentTimeMillis()-startTime)+"ms");
        System.out.println("LevenshteinDistance: " + LevenshteinDistance.Compute(supreme.toString(),DNA.toString()));

    }


}
