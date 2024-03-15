package ccd.algorithms;

//import java.util.BitSet;
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
	 * @param container
	 *            tested if contains other BitSet
	 * @param contained
	 *            tested if contained in other BitSet
	 * @return whether the container BitSet contains the other BitSet
	 */
	public static boolean contains(BitSet container, BitSet contained) {
		BitSet copy = (BitSet) container.clone();
		copy.and(contained);
		return copy.equals(contained);
	}

	/**
	 * Whether the two given BitSets are disjoint.
	 * 
	 * @param first
	 *            tested if disjoint from second
	 * @param second
	 *            tested if disjoint from first
	 * @return whether the two given BitSets are disjoint.
	 */
	public static boolean disjoint(BitSet first, BitSet second) {
		BitSet copy = (BitSet) first.clone();
		copy.and(second);
		return copy.isEmpty();
	}

	/**
	 * Whether the two given BitSets intersect properly, that is, they intersect
	 * and neither is contained in the other.
	 * 
	 * @param first
	 *            tested if intersecting properly with second
	 * @param second
	 *            tested if intersecting properly with first
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

	public static void main(String[] args) {
		BitSet first = new BitSet(5);
		BitSet second = new BitSet(5);
		first.set(0);
		first.set(1);
		second.set(1);
		second.set(2);
		System.out.println("0,1 - 1,2");
		System.out.println("iP: " + intersectProperly(first, second));
		System.out.println("d:  " + disjoint(first, second));
		System.out.println("c:  " + contains(first, second));

		first.set(2);
		System.out.println("0,1,2 - 1,2");
		System.out.println("iP: " + intersectProperly(first, second));
		System.out.println("d:  " + disjoint(first, second));
		System.out.println("c:  " + contains(first, second));
		System.out.println("c-: " + contains(second, first));

		second.clear();
		second.set(3);
		System.out.println("0,1,2 - 3");
		System.out.println("iP: " + intersectProperly(first, second));
		System.out.println("d:  " + disjoint(first, second));
		System.out.println("c:  " + contains(first, second));
		System.out.println("c-: " + contains(second, first));
	}

	public static BitSet getLexicographicFirst(BitSet bitsOne, BitSet bitsTwo) {
		if (bitsOne.equals(bitsTwo)) {
			// if they are equal, then we can return either bitset
			return bitsOne;
		}

		// mask1 will contain all bits (taxa) that are in bitsTwo, but not in bitsOne
		BitSet mask1 = new BitSet();
		mask1.or(bitsOne);
		mask1.andNot(bitsTwo);
		mask1.or(bitsTwo);
		mask1.andNot(bitsOne);

		// mask2 will contain all bits (taxa) that are in bitsOne, but not in bitsTwo
		BitSet mask2 = new BitSet();
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
}
