package ccd.algorithms;

import ccd.model.BitSet;

/**
 * This class provides methods to work with BitSets.
 *
 * @author Jonathan Klawitter
 */
public class BitSetUtil {

    /**
     * Whether the container BitSet contains the other "contained" BitSet.
     *
     * @param container tested if contains other BitSet
     * @param contained tested if contained in other BitSet
     * @return whether the container BitSet contains the other BitSet
     */
    public static boolean contains(BitSet container, BitSet contained) {
        return container.contains(contained);
    }

    /**
     * Whether the two given BitSets are disjoint.
     *
     * @param first  tested if disjoint from second
     * @param second tested if disjoint from first
     * @return whether the two given BitSets are disjoint.
     */
    public static boolean disjoint(BitSet first, BitSet second) {
        return first.disjoint(second);
    }

    /**
     * Whether the two given BitSets intersect properly, that is, they intersect
     * and neither is contained in the other.
     *
     * @param first  tested if intersecting properly with second
     * @param second tested if intersecting properly with first
     * @return whether the two given BitSets intersect properly
     */
    public static boolean intersectProperly(BitSet first, BitSet second) {
        BitSet copy = (BitSet) first.clone();
        copy.and(second);
        if (copy.isEmpty()) {
            return false;
        }
        BitSet copyAnd = (BitSet) copy.clone();
        copy.xor(first);
        copyAnd.xor(second);
        return !copy.isEmpty() && !copyAnd.isEmpty();
    }

    public static BitSet getLexicographicFirst(BitSet bitsOne, BitSet bitsTwo) {
        if (bitsOne.equals(bitsTwo)) {
            // if they are equal, then we can return either bitset
            return bitsOne;
        }

        // mask1 will contain all bits (taxa) that are in bitsTwo, but not in bitsOne
        BitSet mask1 = BitSet.newBitSet(bitsOne.size());
        mask1.or(bitsOne);
        mask1.andNot(bitsTwo);
        mask1.or(bitsTwo);
        mask1.andNot(bitsOne);

        // mask2 will contain all bits (taxa) that are in bitsOne, but not in bitsTwo
        BitSet mask2 = BitSet.newBitSet(bitsOne.size());
        mask2.or(bitsTwo);
        mask2.andNot(bitsOne);
        mask2.or(bitsOne);
        mask2.andNot(bitsTwo);

        // to know which clade is lexico smaller, we only have to compare
        // the position of the first set bits in the masks
        int firstDifferentBitInTwo = mask1.nextSetBit(0);
        int firstDifferentBitInOne = mask2.nextSetBit(0);
        return (firstDifferentBitInOne < firstDifferentBitInTwo) ? bitsOne : bitsTwo;
    }

    /**
     * @param first  bitset
     * @param second bitset
     * @return number of bits set to true in both bitsets
     */
    public static int overlapSize(BitSet first, BitSet second) {
        BitSet copy = (BitSet) first.clone();
        copy.and(second);
        return copy.cardinality();
    }

    /**
     * @param bits to be toggled
     * @param end  highest index + 1 for toggled bits
     * @return BitSet with each bit opposite to what it is in the given bits
     */
    public static BitSet getToggled(BitSet bits, int end) {
        int size = bits.size();
        BitSet mask = BitSet.newBitSet(size);
        mask.set(1, end);
        mask.andNot(bits);

        return mask;
    }
}
