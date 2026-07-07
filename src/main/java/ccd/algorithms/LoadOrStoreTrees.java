package ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.parser.NexusParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

public class LoadOrStoreTrees {

    public static List<Tree> loadTrees(File treeFile, double burninPercentage) throws IOException {
        NexusParser parser = new NexusParser();
        parser.parseFile(treeFile);
        System.out.println("loading tree file: " + treeFile);
        List<Tree> fullTreeList = parser.trees;
        List<Tree> treesToUse;
        if (burninPercentage == 0) {
            treesToUse = fullTreeList;
        } else {
            int numDiscardedTrees = (int) (fullTreeList.size() * burninPercentage);
            int numUsedTrees = fullTreeList.size() - numDiscardedTrees;
            treesToUse = new ArrayList<Tree>(numUsedTrees);
            treesToUse.addAll(fullTreeList.subList(numDiscardedTrees, fullTreeList.size()));
        }
        return treesToUse;
    }

    /**
     * Loads at most {@code sampleCount} trees, thinned evenly across the post-burnin chain, parsing
     * <em>only</em> the sampled trees: the Newick strings of the trees in between are read past in the
     * reader but never handed to the {@link TreeParser}, which is the expensive step. This is the
     * streaming counterpart of {@link #loadTrees(File, double)} followed by even thinning -- it
     * returns the same trees, at the same indices, but skips building the discarded ones.
     *
     * <p>The thinning matches the in-memory idiom: with {@code total} trees in the file, the first
     * {@code floor(total * burnin)} are burn-in; of the remaining {@code usable} trees, tree
     * {@code numDiscarded + floor(i * usable / sampleCount)} is taken for {@code i = 0..sampleCount-1}.
     * If {@code usable <= sampleCount} every post-burnin tree is returned (nothing is thinned).
     *
     * <p>Two streaming passes are made over the file: the first counts the trees (no Newick parsing)
     * to fix {@code total}, the second parses only the selected indices. Both passes use one shared
     * command reader so the indices line up exactly.
     */
    public static List<Tree> loadTrees(File treeFile, double burninPercentage, int sampleCount)
            throws IOException {
        if (sampleCount <= 0) {
            throw new IllegalArgumentException("sampleCount must be positive, got " + sampleCount);
        }
        System.out.println("loading tree file: " + treeFile + " (subsampling up to " + sampleCount
                + " of the post-burnin trees, skipping the rest in the reader)");
        SubsamplingNexusParser parser = new SubsamplingNexusParser();

        // Pass 1: count trees without parsing any Newick, so even thinning has the true total.
        int total = parser.streamTreesBlock(treeFile, null, null);

        int numDiscarded = (int) (total * burninPercentage);
        int usable = total - numDiscarded;
        BitSet wanted = new BitSet(total);
        if (usable <= sampleCount) {
            wanted.set(numDiscarded, total);
        } else {
            double step = (double) usable / sampleCount;
            for (int i = 0; i < sampleCount; i++) {
                wanted.set(numDiscarded + (int) (i * step));
            }
        }

        // Pass 2: parse only the selected trees.
        List<Tree> trees = new ArrayList<>(wanted.cardinality());
        parser.streamTreesBlock(treeFile, wanted, trees);
        return trees;
    }

    /**
     * A {@link NexusParser} that can stream the trees block, parsing only a chosen subset of trees.
     * It reuses the superclass's translate/taxa handling (so a tree set with or without a translate
     * block parses exactly as {@link NexusParser} would) and replicates only the small command reader,
     * which is package-private in {@link NexusParser} and so not otherwise reachable.
     */
    private static final class SubsamplingNexusParser extends NexusParser {

        /**
         * Streams the trees block of {@code treeFile}. Every {@code tree} command is counted; a tree
         * whose zero-based index is set in {@code wanted} is parsed into a {@link Tree} and appended to
         * {@code out}. When {@code wanted} is {@code null} nothing is parsed (a pure count pass).
         *
         * @return the total number of trees in the file
         */
        int streamTreesBlock(File treeFile, BitSet wanted, List<Tree> out) throws IOException {
            try (BufferedReader fin = new BufferedReader(new FileReader(treeFile))) {
                lineNr = 0;
                // Advance to the trees block, parsing any preceding taxa block (needed for taxon names
                // when there is no translate block), exactly as NexusParser.parseFile does.
                String line;
                while ((line = nextLine(fin)) != null) {
                    String lower = line.toLowerCase();
                    if (lower.matches("^\\s*begin\\s+taxa;\\s*$")) {
                        parseTaxaBlock(fin);
                    } else if (lower.matches("^\\s*begin\\s+trees;\\s*$")) {
                        // The header (taxa block, block starts) is tiny; the trees block is the bulk of
                        // a multi-MB file, so read it through a chunked buffer rather than the per-char
                        // (synchronized) BufferedReader.read() that otherwise dominates parse time.
                        return readTrees(new FastCharReader(fin), wanted, out);
                    }
                }
            }
            return 0;
        }

        private int readTrees(FastCharReader fin, BitSet wanted, List<Tree> out) throws IOException {
            int origin = -1;
            List<String> blockTaxa = this.taxa;

            String command = readCommand(fin);
            // optional translate block, parsed once up front
            if (command != null && commandName(command).equals("translate")) {
                Map<String, String> translationMap = parseTranslateCommand(commandArgs(command));
                origin = getIndexedTranslationMapOrigin(translationMap);
                if (origin != -1) {
                    blockTaxa = getIndexedTranslationMap(translationMap, origin);
                }
                command = readCommand(fin);
            }

            int current = 0;
            for (; command != null; command = readCommand(fin)) {
                // Only the command name is needed to classify; the full (whitespace-normalised)
                // arguments are extracted lazily, i.e. only for the trees actually being parsed.
                String name = commandName(command);
                if (name.equals("end")) {
                    break;
                }
                if (name.equals("tree")) {
                    if (wanted != null && wanted.get(current)) {
                        out.add(parseTree(commandArgs(command), blockTaxa, origin, current));
                    }
                    current++;
                }
            }
            return current;
        }

        /** Mirrors NexusParser.parseTreesBlock's per-tree handling for a single selected tree. */
        private Tree parseTree(String treeString, List<String> blockTaxa, int origin, int current) {
            int i = treeString.indexOf('(');
            String id = "" + current;
            try {
                id = treeString.substring(5, i).split("=")[0].trim();
            } catch (Exception e) {
                // ignore, keep the index as id
            }
            if (i > 0) {
                treeString = treeString.substring(i);
            }
            TreeParser treeParser;
            if (origin != -1) {
                treeParser = new TreeParser(blockTaxa, treeString, origin, false);
            } else {
                try {
                    treeParser = new TreeParser(blockTaxa, treeString, 0, false);
                } catch (ArrayIndexOutOfBoundsException e) {
                    treeParser = new TreeParser(blockTaxa, treeString, 1, false);
                }
            }
            treeParser.setID(id);
            return treeParser;
        }

        /* ---- replicated command reader (NexusParser's is package-private) ---- */

        /**
         * Reads the next nexus command: characters up to the terminating {@code ;}, with nexus
         * comments ({@code [...]}) and quoted strings read through verbatim. Returns {@code null} at
         * end of file. Mirrors {@link NexusParser}'s {@code readNextCommand}.
         */
        private String readCommand(FastCharReader fin) throws IOException {
            StringBuilder sb = new StringBuilder();
            int nextVal;
            while ((nextVal = fin.read()) >= 0) {
                char c = (char) nextVal;
                if (c == ';') {
                    break;
                }
                sb.append(c);
                switch (c) {
                    case '[' -> readComment(fin, sb);
                    case '"', '\'' -> readString(fin, sb, c);
                    case '\n' -> lineNr += 1;
                    default -> {
                    }
                }
            }
            return sb.toString().isEmpty() ? null : sb.toString();
        }

        private void readComment(FastCharReader fin, StringBuilder sb) throws IOException {
            int nextVal;
            while ((nextVal = fin.read()) >= 0) {
                char c = (char) nextVal;
                sb.append(c);
                if (c == ']') {
                    return;
                }
                if (c == '"' || c == '\'') {
                    readString(fin, sb, c);
                }
                if (c == '\n') {
                    lineNr += 1;
                }
            }
            throw new IOException("Unterminated comment.");
        }

        private void readString(FastCharReader fin, StringBuilder sb, char delim) throws IOException {
            int nextVal;
            while ((nextVal = fin.read()) >= 0) {
                char c = (char) nextVal;
                sb.append(c);
                if (c == delim) {
                    return;
                }
                if (c == '\n') {
                    lineNr += 1;
                }
            }
            throw new IOException("Unterminated string.");
        }

        /**
         * Pulls characters from an underlying {@link Reader} through a large in-memory buffer. The
         * nexus command reader consumes the trees block one character at a time; doing that straight
         * off a {@link BufferedReader} pays the per-call cost of its {@code synchronized} {@code read()}
         * (plus bounds work) on every character, which dominated parse time on multi-MB tree files.
         * Here one {@code read(char[])} refills ~64k characters, so that cost is paid ~once per 64k
         * characters instead of once per character. Behaviour is otherwise identical: same characters,
         * same order. Not thread-safe (one instance per parse pass).
         */
        private static final class FastCharReader {
            private final Reader in;
            private final char[] buf = new char[1 << 16];
            private int pos = 0;
            private int len = 0;

            FastCharReader(Reader in) {
                this.in = in;
            }

            /** @return the next character, or -1 at end of stream. */
            int read() throws IOException {
                if (pos >= len) {
                    len = in.read(buf);
                    pos = 0;
                    if (len <= 0) {
                        return -1;
                    }
                }
                return buf[pos++];
            }
        }

        /**
         * Command name: the first whitespace-delimited token, lower-cased -- matching the {@code command}
         * field of {@link NexusParser}'s NexusCommand. Scanned directly (no whitespace normalisation of
         * the whole command) so classifying a skipped tree does not touch its Newick string.
         */
        private static String commandName(String command) {
            int n = command.length(), i = 0;
            while (i < n && Character.isWhitespace(command.charAt(i))) i++;
            int start = i;
            while (i < n && !Character.isWhitespace(command.charAt(i))) i++;
            return command.substring(start, i).toLowerCase();
        }

        /**
         * Command arguments: everything after the command token, with internal whitespace collapsed to
         * single spaces -- matching the {@code arguments} field of {@link NexusParser}'s NexusCommand.
         */
        private static String commandArgs(String command) {
            String trimmed = command.trim().replaceAll("\\s+", " ");
            String name = trimmed.split(" ")[0];
            return trimmed.length() > name.length() + 1 ? trimmed.substring(name.length() + 1) : "";
        }
    }

    // TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.getAbsolutePath(), 10);
    // AbstractCCD initialCCD1 = new CCD1(treeSet, false);  // Create initial CCD1 object

    public static void storeTrees(List<Tree> trees, File treeFile) throws FileNotFoundException {
        try (PrintStream handle = new PrintStream(treeFile)) {

            trees.get(0).init(handle); // Write header info
            handle.println(); // Ensure the semicolon before the first tree is on its own line

            for (int i = 0; i < trees.size(); i++) {
                Tree tree = trees.get(i);

                int id = (tree.getID() == null) ? i : Integer.parseInt(tree.getID().replace("_", ""));
                tree.log(id, handle);

                handle.println(); // blank line between trees
            }

            trees.get(trees.size() - 1).close(handle);
        }
    }
}
