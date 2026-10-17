package ccd.model.bitsets;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Copies of a BitSet must keep the source's capacity, not shrink to its highest set bit.
 *
 * <p>The bitwise operations iterate over {@code this.words.length} and index the operand directly,
 * so an undersized copy makes them throw. This bites only above 256 bits, where the generic BitSet
 * is used instead of the fixed-size subclasses, and only when the top words are empty -- the common
 * case for a clade that does not contain the highest-numbered taxa.
 */
public class BitSetCopyTest {

    @Test
    public void copyKeepsCapacityAboveTheSpecialisedSizes() {
        for (int nbits : new int[]{276, 320, 512, 1000}) {
            BitSet full = BitSet.newBitSet(nbits);
            full.set(nbits - 1);

            BitSet sparse = BitSet.newBitSet(nbits);
            sparse.set(3);
            BitSet copy = BitSet.newBitSet(sparse);

            assertEquals(sparse.size(), copy.size(),
                    "copy must keep the source capacity at " + nbits + " bits");

            BitSet a = BitSet.newBitSet(full);
            a.andNot(copy);
            assertTrue(a.get(nbits - 1), "andNot must not clear unrelated high bits");

            BitSet b = BitSet.newBitSet(copy);
            b.andNot(full);
            assertTrue(b.get(3), "andNot must keep the low bit");
            assertFalse(copy.intersects(full), "disjoint sets must not intersect");
        }
    }

    @Test
    public void copyOfAnEmptySetIsUsable() {
        BitSet empty = BitSet.newBitSet(400);
        BitSet copy = BitSet.newBitSet(empty);
        assertEquals(empty.size(), copy.size(), "an empty copy must still have capacity");
        BitSet other = BitSet.newBitSet(400);
        other.set(399);
        copy.andNot(other);
        copy.or(other);
        assertTrue(copy.get(399));
    }
}
