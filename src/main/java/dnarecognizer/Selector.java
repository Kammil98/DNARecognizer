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
            for(selected = 0; randomVal >= 0; selected++){
                randomVal -= population.get(selected).getFitVal();
            }
            selected--; // one down to choose correct DNA chain
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
        ArrayList<DNAChain> indexes = new ArrayList<>();
        ArrayList<DNAChain> newPopulation = new ArrayList<>();
        DNAChain winner;
        int currPopSize, contestGroupSize = 4;
        currPopSize = population.size();
        while(newPopulation.size() < DEMAND_SIZE){
            //choosing indexes of members of population to contest
            for(int i = 0; i < contestGroupSize ; i++) {
                indexes.add(population.remove(DNARecognizer.getGENERATOR().nextInt(currPopSize)));
                currPopSize--;
            }

            //choosing index of best in group
            winner = indexes.get(0);
            for(int i = 1; i < contestGroupSize ; i++){
                if(indexes.get(i).getFitVal() > winner.getFitVal()){
                    winner = indexes.get(i);
                }
            }
            newPopulation.add(winner);
            indexes.clear();
        }
        population.addAll(newPopulation);
    }

    /**
     * cut population to demandSize best
     * DNA Chains
     * @param population group of DNAChains, which will take part in selection
     */
    public void ranking(ArrayList<DNAChain> population){ //probably not working | sort what? size/fitValue?
        population.sort(Collections.reverseOrder());
        population.subList(0, DEMAND_SIZE);
    }
}
