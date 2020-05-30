package dnarecognizer;

import dnarecognizer.models.DNAChain;
import dnarecognizer.models.Oligonucleotide;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Collections;

import static org.junit.jupiter.api.Assertions.*;

class DNARecognizerTest {

    public DNAChain member1, member2, oldMember1, oldMember2;

    @BeforeEach
    void init() {
        member1 = new DNAChain();
        member2 = new DNAChain();
        oldMember1 = new DNAChain();
        oldMember2 = new DNAChain();
        Oligonucleotide olig;
        for(int i = 0; i < 10; i++){
            olig = new Oligonucleotide(new StringBuilder("AACCTTGG"));
            member1.getOligonucleotides().add(olig);
        }
        member2.setOligonucleotides(new ArrayList<>(member1.getOligonucleotides()));
        Collections.shuffle(member2.getOligonucleotides());
        oldMember1.setOligonucleotides(new ArrayList<>(member1.getOligonucleotides()));
        oldMember2.setOligonucleotides(new ArrayList<>(member2.getOligonucleotides()));
    }

    @Test
    void mutateMemberByShuffling() {
        System.out.println("mutateMemberByShuffling");
        double mutationSize = 20.0d;
        DNARecognizer.mutateMemberByShuffling(member1, mutationSize);
        assertNotEquals(oldMember1, member1, "mutateMemberByShuffling did not " +
                "change order. Can be accident(probability) try to run again");
    }

    @Test
    void mutateMemberByAdding() {
        System.out.println("mutateMemberByAdding: check if oligs were added");
        double mutationSize = 20.0d;
        int amount = (int) Math.ceil(mutationSize *
                (member1.getOligonucleotides().size() / 100.0d));
        DNARecognizer.mutateMemberByAdding(member1, mutationSize);
        member1.getOligonucleotides().removeAll(oldMember1.getOligonucleotides());
        assertEquals(member1.getOligonucleotides().size(), amount, "mutateMemberByAdding did not add Oligonucleotide.");
    }

    @Test
    void crossTwoMembers() {
        System.out.println("crossTwoMembers : check is all oligs copied");
        DNAChain child = DNARecognizer.crossTwoMembers(member1, member2, 50.0d);
        for(Oligonucleotide olig: member1.getOligonucleotides()){
            if(!child.getOligonucleotides().contains(olig))
                fail("crossTwoMembers did not copied All oligs from mother.");
        }
    }

}