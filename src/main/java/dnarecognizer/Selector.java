package dnarecognizer;

import dnarecognizer.models.DNAChain;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author kamil_2
 */
public class Selector {

    //size, to which all functions will cut population
    private final int DEMAND_SIZE;
    public Selector(int demandSize){
        this.DEMAND_SIZE = demandSize;
    }

    /**
     * select best of group
     * @param group group of DNAChains, which will take part in selection
     */
    public DNAChain bestOfGroup(ArrayList<DNAChain> group){
        int winner = 0;
        for(int i = 1; i < group.size() ; i++)
            if(group.get(i).getFitVal() > group.get(winner).getFitVal())
                winner = i;
        return group.get(winner);
    }

    /**
     * cut population to demandSize size
     * remove random elements with probability
     * based on fitValue
     * @param population group of DNAChains, which will take part in selection
     */
    public void roulette(ArrayList<DNAChain> population){
        int sum = 0, selected, randomVal;

        //stream calculate sum of all fitValues
        sum = population.stream()
                .map(DNAChain::getFitVal)
                .reduce(sum, Integer::sum);

        while(population.size() > DEMAND_SIZE){
            randomVal = (int) (sum * DNARecognizer.getGENERATOR().nextDouble());
            for(selected = 0; randomVal > 0; selected++){
                randomVal -= population.get(selected).getFitVal();
            }
            //removing element and update sum
            sum -= population.get(selected).getFitVal();
            population.remove(selected);
        }
    }

    /**
     * cut population to demandSize size
     * take group of random members from population
     * and remove all above best in group
     * @param population group of DNAChains, which will take part in selection
     */
    public void contest(ArrayList<DNAChain> population){
        int[] indexes = new int[4];
        int currPopSize, contestGroupSize = 4;
        int winner, looser;
        while(population.size() > DEMAND_SIZE){
            //choosing indexes of members of population to contest
            currPopSize = population.size();
            for(int i = 0; i < contestGroupSize ; i++)
                indexes[i] = DNARecognizer.getGENERATOR().nextInt(currPopSize);

            //choosing index of best in group
            winner = indexes[0];
            for(int i = 1; i < contestGroupSize ; i++){
                if(population.get(i).getFitVal() > population.get(winner).getFitVal()){
                    looser = winner;
                    winner = i;
                }
                else
                    looser = i;
                if(population.size() > DEMAND_SIZE)
                    population.remove(looser);
            }
        }
    }

    /**
     * cut population to demandSize best
     * DNA Chains
     * @param population group of DNAChains, which will take part in selection
     */
    public void ranking(ArrayList<DNAChain> population){
        Collections.sort(population);
        population.subList(0, DEMAND_SIZE);
    }
}
