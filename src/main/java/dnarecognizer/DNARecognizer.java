package dnarecognizer;

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

    //in range 300-1000 nucleotides
    private static final int DNA_SIZE = 300;
    // in range 8-10 nucleotides
    private static final int OLIGONUCLEOTIDE_SIZE = 8;
    private static final int POPULATION_SIZE = 100;
    private static final int EXPANDED_POPULATION_SIZE = 150;
    private static final int CHILDREN_NO_PER_ONE_MATCH = 4;
    private static final double COPIED_PERCENT_OF_DNA = 40.0d;
    private static final Random GENERATOR = new Random();
    private static final Selector SELECTOR = new Selector(POPULATION_SIZE);

    /**
     * selecting and crossing members of population,
     * then adding results to our population
     * @param population all possible results of our metaheuristic
     */
    public static void crossover(ArrayList<DNAChain> population) {
        int breedNo = EXPANDED_POPULATION_SIZE - population.size();
        int groupSize = POPULATION_SIZE / 20;
        ArrayList<DNAChain> group = new ArrayList<>();
        DNAChain member1, member2;
        while(breedNo > 0){

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
            if(CHILDREN_NO_PER_ONE_MATCH % 2 == 0)
                population.add(crossTwoMembers(member1, member2, COPIED_PERCENT_OF_DNA));

            breedNo -= CHILDREN_NO_PER_ONE_MATCH;
        }
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

        fatherPartBegin = GENERATOR.nextInt(father.getOligonucleotides().size() - copiedPartAmount);
        fatherPartEnd = fatherPartBegin + copiedPartAmount;
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
        if (optimumNo < 3)
            return;
        Iterator<DNAChain> iter = population.iterator();
        for(DNAChain member = iter.next();
            iter.hasNext() && optimumNo > 3;
            member = iter.next()){
            if(optimumNo > POPULATION_SIZE / 5)
                shufflingDegree = Math.max(100.0d, shufflingDegree + 70.0d);
            else
                shufflingDegree = Math.max(100.0d, shufflingDegree + 30.0d);
            mutateMemberByShuffling(member, shufflingDegree);
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
                oliginucleotides.size() - shufflingAmount);
        List<Oligonucleotide> oligToShuffle = oliginucleotides
                .subList(shufflingBegin, shufflingBegin + shufflingAmount);
        Collections.shuffle(oligToShuffle);
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
        if(progress < 50){
            SELECTOR.roulette(population);
        }
        else if(progress < 90){
            SELECTOR.contest(population);
        }
        else{
            SELECTOR.ranking(population);
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        ArrayList<DNAChain> population = new ArrayList<>(EXPANDED_POPULATION_SIZE);

    }



}
