package ccd.model;


public class BitSet extends java.util.BitSet {

//    /**
//     * The internal field corresponding to the serialField "bits".
//     */
//    private long[] words;
//
//    /**
//     * The number of words in the logical size of this BitSet.
//     */
//    private transient int wordsInUse = 0;
//
//    /**
//     * Whether the size of "words" is user-specified.  If so, we assume
//     * the user knows what he's doing and try harder to preserve it.
//     */
//    private transient boolean sizeIsSticky = false;

	protected BitSet() {
		super();
	}
	
	/** not to be used for production code **/
	@Deprecated
	public static BitSet newBitSetForTesting(int nbits) {
		return new BitSet(nbits);
	}
	public static BitSet newBitSet(int nbits) {
		// if (true) return new BitSet(nbits);
		if (nbits <= 64) {
			return new BitSet64();
		} else if (nbits <= 128) {
			return new BitSet128();
		} else if (nbits <= 192) {
			return new BitSet192();
		} else if (nbits <= 256) {
			return new BitSet256();
		} else { 
			return new BitSet(nbits);
		}
	}

	public static BitSet newBitSet(BitSet other) {
		if (other instanceof BitSet128 set) {
			return new BitSet128(set);
		}
		if (other instanceof BitSet192 set) {
			return new BitSet192(set);
		}
		if (other instanceof BitSet64 set) {
			return new BitSet64(set);
		}
		if (other instanceof BitSet256 set) {
			return new BitSet256(set);
		}
		BitSet b = new BitSet(other.length());
		b.or(other);
		return b;
	}
	

	private BitSet(int nbits) {
		super(nbits);
//        // nbits can't be negative; size 0 is OK
//        if (nbits < 0)
//            throw new NegativeArraySizeException("nbits < 0: " + nbits);
//
//        words = new long[wordIndex(nbits-1) + 1];
//        sizeIsSticky = true;
	}
	
//	   /**
//     * Sets the bit at the specified index to {@code true}.
//     *
//     * @param  bitIndex a bit index
//     * @throws IndexOutOfBoundsException if the specified index is negative
//     * @since  1.0
//     */
    public void set(int bitIndex) {
//		throw new RuntimeException();
    	super.set(bitIndex);
//        if (bitIndex < 0)
//            throw new IndexOutOfBoundsException("bitIndex < 0: " + bitIndex);
//
//        int wordIndex = wordIndex(bitIndex);
//        expandTo(wordIndex);
//
//        words[wordIndex] |= (1L << bitIndex); // Restores invariants
//
//        checkInvariants();
    }
//
//    /**
//     * Sets the bits from the specified {@code fromIndex} (inclusive) to the
//     * specified {@code toIndex} (exclusive) to {@code true}.
//     *
//     * @param  fromIndex index of the first bit to be set
//     * @param  toIndex index after the last bit to be set
//     * @throws IndexOutOfBoundsException if {@code fromIndex} is negative,
//     *         or {@code toIndex} is negative, or {@code fromIndex} is
//     *         larger than {@code toIndex}
//     * @since  1.4
//     */
	public void set(int fromIndex, int toIndex) {
//		throw new RuntimeException();
		super.set(fromIndex, toIndex);
//        checkRange(fromIndex, toIndex);
//
//        if (fromIndex == toIndex)
//            return;
//
//        // Increase capacity if necessary
//        int startWordIndex = wordIndex(fromIndex);
//        int endWordIndex   = wordIndex(toIndex - 1);
//        expandTo(endWordIndex);
//
//        long firstWordMask = WORD_MASK << fromIndex;
//        long lastWordMask  = WORD_MASK >>> -toIndex;
//        if (startWordIndex == endWordIndex) {
//            // Case 1: One word
//            words[startWordIndex] |= (firstWordMask & lastWordMask);
//        } else {
//            // Case 2: Multiple words
//            // Handle first word
//            words[startWordIndex] |= firstWordMask;
//
//            // Handle intermediate words, if any
//            for (int i = startWordIndex+1; i < endWordIndex; i++)
//                words[i] = WORD_MASK;
//
//            // Handle last word (restores invariants)
//            words[endWordIndex] |= lastWordMask;
//        }
//
//        checkInvariants();
	}
//
//    /**
//     * Performs a logical <b>OR</b> of this bit set with the bit set
//     * argument. This bit set is modified so that a bit in it has the
//     * value {@code true} if and only if it either already had the
//     * value {@code true} or the corresponding bit in the bit set
//     * argument has the value {@code true}.
//     *
//     * @param set a bit set
//     */
    public void or(BitSet set) {
//		throw new RuntimeException();
    	super.or(set);
//        if (this == set)
//            return;
//
//        int wordsInCommon = Math.min(wordsInUse, set.wordsInUse);
//
//        if (wordsInUse < set.wordsInUse) {
//            ensureCapacity(set.wordsInUse);
//            wordsInUse = set.wordsInUse;
//        }
//
//        // Perform logical OR on words in common
//        for (int i = 0; i < wordsInCommon; i++)
//            words[i] |= set.words[i];
//
//        // Copy any remaining words
//        if (wordsInCommon < set.wordsInUse)
//            System.arraycopy(set.words, wordsInCommon,
//                             words, wordsInCommon,
//                             wordsInUse - wordsInCommon);
//
//        // recalculateWordsInUse() is unnecessary
//        checkInvariants();
    }
//
//    /**
//     * Returns the index of the first bit that is set to {@code true}
//     * that occurs on or after the specified starting index. If no such
//     * bit exists then {@code -1} is returned.
//     *
//     * <p>To iterate over the {@code true} bits in a {@code BitSet},
//     * use the following loop:
//     *
//     *  <pre> {@code
//     * for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
//     *     // operate on index i here
//     *     if (i == Integer.MAX_VALUE) {
//     *         break; // or (i+1) would overflow
//     *     }
//     * }}</pre>
//     *
//     * @param  fromIndex the index to start checking from (inclusive)
//     * @return the index of the next set bit, or {@code -1} if there
//     *         is no such bit
//     * @throws IndexOutOfBoundsException if the specified index is negative
//     * @since  1.4
//     */
    public int nextSetBit(int fromIndex) {
//		throw new RuntimeException();
    	return super.nextSetBit(fromIndex);
//        if (fromIndex < 0)
//            throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);
//
//        checkInvariants();
//
//        int u = wordIndex(fromIndex);
//        if (u >= wordsInUse)
//            return -1;
//
//        long word = words[u] & (WORD_MASK << fromIndex);
//
//        while (true) {
//            if (word != 0)
//                return (u * BITS_PER_WORD) + Long.numberOfTrailingZeros(word);
//            if (++u == wordsInUse)
//                return -1;
//            word = words[u];
//        }
    }
    /**
     * Returns the index of the first bit that is set to {@code false}
     * that occurs on or after the specified starting index.
     *
     * @param  fromIndex the index to start checking from (inclusive)
     * @return the index of the next clear bit
     * @throws IndexOutOfBoundsException if the specified index is negative
     * @since  1.4
     */
    public int nextClearBit(int fromIndex) {
    	return super.nextClearBit(fromIndex);
//		throw new RuntimeException();
    }
    	
//    /**
//     * Cloning this {@code BitSet} produces a new {@code BitSet}
//     * that is equal to it.
//     * The clone of the bit set is another bit set that has exactly the
//     * same bits set to {@code true} as this bit set.
//     *
//     * @return a clone of this bit set
//     * @see    #size()
//     */
    public Object clone() {
//		throw new RuntimeException();
    	return super.clone();
//        if (! sizeIsSticky)
//            trimToSize();
//
//        try {
//            BitSet result = (BitSet) super.clone();
//            result.words = words.clone();
//            result.checkInvariants();
//            return result;
//        } catch (CloneNotSupportedException e) {
//            throw new InternalError(e);
//        }
    }
    
//    /**
//     * Performs a logical <b>AND</b> of this target bit set with the
//     * argument bit set. This bit set is modified so that each bit in it
//     * has the value {@code true} if and only if it both initially
//     * had the value {@code true} and the corresponding bit in the
//     * bit set argument also had the value {@code true}.
//     *
//     * @param set a bit set
//     */
    public void and(BitSet set) {
//		throw new RuntimeException();
    	super.and(set);
//        if (this == set)
//            return;
//
//        while (wordsInUse > set.wordsInUse)
//            words[--wordsInUse] = 0;
//
//        // Perform logical AND on words in common
//        for (int i = 0; i < wordsInUse; i++)
//            words[i] &= set.words[i];
//
//        recalculateWordsInUse();
//        checkInvariants();
    }
//    
//    /**
//     * Clears all of the bits in this {@code BitSet} whose corresponding
//     * bit is set in the specified {@code BitSet}.
//     *
//     * @param  set the {@code BitSet} with which to mask this
//     *         {@code BitSet}
//     * @since  1.2
//     */
    public void andNot(BitSet set) {
//		throw new RuntimeException();
    	super.andNot(set);
//        // Perform logical (a & !b) on words in common
//        for (int i = Math.min(wordsInUse, set.wordsInUse) - 1; i >= 0; i--)
//            words[i] &= ~set.words[i];
//
//        recalculateWordsInUse();
//        checkInvariants();
    }
//    
//    /**
//     * Performs a logical <b>XOR</b> of this bit set with the bit set
//     * argument. This bit set is modified so that a bit in it has the
//     * value {@code true} if and only if one of the following
//     * statements holds:
//     * <ul>
//     * <li>The bit initially has the value {@code true}, and the
//     *     corresponding bit in the argument has the value {@code false}.
//     * <li>The bit initially has the value {@code false}, and the
//     *     corresponding bit in the argument has the value {@code true}.
//     * </ul>
//     *
//     * @param  set a bit set
//     */
    public void xor(BitSet set) {
//		throw new RuntimeException();
    	super.xor(set);
//        int wordsInCommon = Math.min(wordsInUse, set.wordsInUse);
//
//        if (wordsInUse < set.wordsInUse) {
//            ensureCapacity(set.wordsInUse);
//            wordsInUse = set.wordsInUse;
//        }
//
//        // Perform logical XOR on words in common
//        for (int i = 0; i < wordsInCommon; i++)
//            words[i] ^= set.words[i];
//
//        // Copy any remaining words
//        if (wordsInCommon < set.wordsInUse)
//            System.arraycopy(set.words, wordsInCommon,
//                             words, wordsInCommon,
//                             set.wordsInUse - wordsInCommon);
//
//        recalculateWordsInUse();
//        checkInvariants();
    }
//    
//    /**
//     * Returns true if the specified {@code BitSet} has any bits set to
//     * {@code true} that are also set to {@code true} in this {@code BitSet}.
//     *
//     * @param  set {@code BitSet} to intersect with
//     * @return boolean indicating whether this {@code BitSet} intersects
//     *         the specified {@code BitSet}
//     * @since  1.4
//     */
    public boolean intersects(BitSet set) {
//		throw new RuntimeException();
    	return super.intersects(set);
//        for (int i = Math.min(wordsInUse, set.wordsInUse) - 1; i >= 0; i--)
//            if ((words[i] & set.words[i]) != 0)
//                return true;
//        return false;
    }
//    
//    /**
//     * Returns true if this {@code BitSet} contains no bits that are set
//     * to {@code true}.
//     *
//     * @return boolean indicating whether this {@code BitSet} is empty
//     * @since  1.4
//     */
    public boolean isEmpty() {
//		throw new RuntimeException();
    	return super.isEmpty();
//        return wordsInUse == 0;
    }
//    
//    /**
//     * Returns the number of bits set to {@code true} in this {@code BitSet}.
//     *
//     * @return the number of bits set to {@code true} in this {@code BitSet}
//     * @since  1.4
//     */
    public int cardinality() {
//		throw new RuntimeException();
    	return super.cardinality();
//        int sum = 0;
//        for (int i = 0; i < wordsInUse; i++)
//            sum += Long.bitCount(words[i]);
//        return sum;
    }
//    
//    /**
//     * Returns the number of bits of space actually in use by this
//     * {@code BitSet} to represent bit values.
//     * The maximum element in the set is the size - 1st element.
//     *
//     * @return the number of bits currently in this bit set
//     */
    public int size() {
//		throw new RuntimeException();
    	return super.size();
//        return words.length * BITS_PER_WORD;
    }
//    
//    /**
//     * Returns the "logical size" of this {@code BitSet}: the index of
//     * the highest set bit in the {@code BitSet} plus one. Returns zero
//     * if the {@code BitSet} contains no set bits.
//     *
//     * @return the logical size of this {@code BitSet}
//     * @since  1.2
//     */
    public int length() {
//		throw new RuntimeException();
    	return super.length();
//        if (wordsInUse == 0)
//            return 0;
//
//        return BITS_PER_WORD * (wordsInUse - 1) +
//            (BITS_PER_WORD - Long.numberOfLeadingZeros(words[wordsInUse - 1]));
    }
//    
//    /**
//     * Sets all of the bits in this BitSet to {@code false}.
//     *
//     * @since 1.4
//     */
    public void clear() {
//		throw new RuntimeException();
    	super.clear();
//        while (wordsInUse > 0)
//            words[--wordsInUse] = 0;
    }
//
//
//
//
//
//
//
//    /*
//     * BitSets are packed into arrays of "words."  Currently a word is
//     * a long, which consists of 64 bits, requiring 6 address bits.
//     * The choice of word size is determined purely by performance concerns.
//     */
//    protected static final int ADDRESS_BITS_PER_WORD = 6;
//    protected static final int BITS_PER_WORD = 1 << ADDRESS_BITS_PER_WORD;
//    protected static final int BIT_INDEX_MASK = BITS_PER_WORD - 1;
//
//    /* Used to shift left or right for a partial word mask */
//    protected static final long WORD_MASK = 0xffffffffffffffffL;
//    /**
//     * Given a bit index, return word index containing it.
//     */
//    private static int wordIndex(int bitIndex) {
//        return bitIndex >> ADDRESS_BITS_PER_WORD;
//    }
//
//    private void expandTo(int wordIndex) {
//        int wordsRequired = wordIndex+1;
//        if (wordsInUse < wordsRequired) {
//            ensureCapacity(wordsRequired);
//            wordsInUse = wordsRequired;
//        }
//    }
//    private void ensureCapacity(int wordsRequired) {
//        if (words.length < wordsRequired) {
//            // Allocate larger of doubled size or required size
//            int request = Math.max(2 * words.length, wordsRequired);
//            words = Arrays.copyOf(words, request);
//            sizeIsSticky = false;
//        }
//    }
//
//    /**
//     * Every public method must preserve these invariants.
//     */
//    private void checkInvariants() {
//        assert(wordsInUse == 0 || words[wordsInUse - 1] != 0);
//        assert(wordsInUse >= 0 && wordsInUse <= words.length);
//        assert(wordsInUse == words.length || words[wordsInUse] == 0);
//    }
//
//    private static void checkRange(int fromIndex, int toIndex) {
//        if (fromIndex < 0)
//            throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);
//        if (toIndex < 0)
//            throw new IndexOutOfBoundsException("toIndex < 0: " + toIndex);
//        if (fromIndex > toIndex)
//            throw new IndexOutOfBoundsException("fromIndex: " + fromIndex +
//                                                " > toIndex: " + toIndex);
//    }
//    private void recalculateWordsInUse() {
//        // Traverse the bitset until a used word is found
//        int i;
//        for (i = wordsInUse-1; i >= 0; i--)
//            if (words[i] != 0)
//                break;
//
//        wordsInUse = i+1; // The new logical size
//    }
//    private void trimToSize() {
//        if (wordsInUse != words.length) {
//            words = Arrays.copyOf(words, wordsInUse);
//            checkInvariants();
//        }
//    }
//
//    public int hashCode() {
//        long h = 1234;
//        for (int i = wordsInUse; --i >= 0; )
//            h ^= words[i] * (i + 1);
//
//        return (int)((h >> 32) ^ h);
//    }
    
    
    public String toString() {

        final int MAX_INITIAL_CAPACITY = Integer.MAX_VALUE - 8;
        int numBits =  cardinality();
        // Avoid overflow in the case of a humongous numBits
        int initialCapacity = (numBits <= (MAX_INITIAL_CAPACITY - 2) / 6) ?
            6 * numBits + 2 : MAX_INITIAL_CAPACITY;
        StringBuilder b = new StringBuilder(initialCapacity);
        b.append('{');

        int i = nextSetBit(0);
        if (i != -1) {
            b.append(i);
            while (true) {
                if (++i < 0) break;
                if ((i = nextSetBit(i)) < 0) break;
                int endOfRun = nextClearBit(i);
                do { b.append(", ").append(i); }
                while (++i != endOfRun);
            }
        }

        b.append('}');
        return b.toString();
    }

    /* returns whether this BitSet contains the other **/
    public boolean contains(BitSet other) {
		BitSet copy = (BitSet) this.clone();
		copy.and(other);
		return copy.equals(other);
    }

	/**returns whether this BitSet is disjoint from the other **/
    public boolean disjoint(BitSet other) {
		BitSet copy = (BitSet) this.clone();
		copy.and(other);
		return copy.isEmpty();
    }
}
