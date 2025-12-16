package ccd.algorithms.regularisation;

import ccd.model.AbstractCCD;
import ccd.model.Clade;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

public class CCDExpansion {

    int maxExpansionFactor;

    public CCDExpansion() {
        this(-1);
    }

    public CCDExpansion(int maxExpansionFactor) {
        this.maxExpansionFactor = maxExpansionFactor;
    }

    public void expandCCD(AbstractCCD ccd) {
        if (maxExpansionFactor == 0) return;

        Stream<Clade> cladesToExpand = ccd.getClades().stream();

        // 1. we only expand the most frequent clades if necessary
        if (this.maxExpansionFactor != -1) {
            int numCladesToConsider = maxExpansionFactor * ccd.getNumberOfLeaves();
            cladesToExpand = cladesToExpand.sorted(
                    Comparator.comparingInt(x -> -x.getNumberOfOccurrences())
            ).limit(numCladesToConsider);
        }

        // 2. sort clades
        List<Clade> clades = cladesToExpand.sorted(Comparator.comparingInt(x -> x.size())).toList();

        // 2. clade buckets
        // for easier matching of child clades, we want to group them by size

        // 2.i init clade buckets
        List<Set<Clade>> cladeBuckets = new ArrayList<Set<Clade>>(ccd.getNumberOfLeaves());
        for (int i = 0; i < ccd.getNumberOfLeaves(); i++) {
            cladeBuckets.add(new HashSet<Clade>());
        }

        // 2.ii fill clade buckets
        for (Clade clade : clades) {
            cladeBuckets.get(clade.size() - 1).add(clade);
        }

        // 3. find missing clade partitions
        this.findChildPartitions(clades, cladeBuckets);
    }

    private void findChildPartitions(List<Clade> parentClades, List<Set<Clade>> cladeBuckets) {
        BitSet helperBits = BitSet.newBitSet(parentClades.get(0).getCCD().getSizeOfLeavesArray());
        int total = parentClades.size();
        int count = 0;

        for (Clade parent : parentClades) {
            count++;
            if (count % 2000 == 0) {
                System.out.println("Progress: " + count + " / " + total + " parent clades processed");
            }
            this.findChildPartitionsOf(parent, helperBits, cladeBuckets);
        }
    }

    private void findChildPartitionsOf(Clade parent, BitSet helperBits, List<Set<Clade>> cladeBuckets) {
        // we skip leaves and cherries as they have no/only one partition
        if (parent.isLeaf() || parent.isCherry()) {
            return;
        }

        BitSet parentBits = parent.getCladeInBits();

        // otherwise we check if we find a larger partner clade for any
        // smaller clade that together partition the parent clade;
        for (int j = 1; j <= parent.size() / 2; j++) {
            for (Clade smallChild : cladeBuckets.get(j - 1)) {
                BitSet smallChildBits = smallChild.getCladeInBits();
                this.findPartitionHelper(smallChild, parent, helperBits, parentBits, smallChildBits);
            }
        }
    }

    private void findPartitionHelper(Clade child, Clade parent, BitSet helperBits, BitSet parentBits, BitSet childBits) {
        // check whether child clade is contained in parent clade
        helperBits.clear();
        helperBits.or(parentBits);
        helperBits.and(childBits);

        if (helperBits.equals(childBits) && !child.getParentClades().contains(parent)) {
            // here helperBits equal childBits, so with an XOR
            // with the parentBits we get the bits of the potential partner clade
            helperBits.xor(parentBits);
            Clade otherChild = child.getCCD().getCladeMapping().get(helperBits);
            if (otherChild != null) {
                parent.createCladePartition(child, otherChild);
            }
        }
    }

    @Override
    public String toString() {
        return switch (this.maxExpansionFactor) {
            case (-1) -> "full expansion";
            case (0) -> "no expansion";
            default -> "expansion with factor " + this.maxExpansionFactor;
        };
    }
}
