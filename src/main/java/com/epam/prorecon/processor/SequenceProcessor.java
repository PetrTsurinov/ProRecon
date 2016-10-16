package com.epam.prorecon.processor;

import com.epam.prorecon.entity.ExonPosition;
import com.epam.prorecon.entity.ExonType;
import htsjdk.variant.variantcontext.VariantContext;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

public class SequenceProcessor {
    private static final Logger LOGGER = LoggerFactory.getLogger(SequenceProcessor.class);

    private static final Map<Character, Character> COMPLEMENT_MRNA_NUCLEOTIDE_MAP =
            Collections.unmodifiableMap(new HashMap<Character, Character>() {{
                put('A', 'U'); put('T', 'A'); put('C', 'G'); put('G', 'C');
            }});
    private static final Map<String, String> NUCLEOTIDE_TRIPLET_AMINO_ACID_MAP =
            Collections.unmodifiableMap(new HashMap<String, String>() {{
                put("UUU", "F"); put("UUC", "F"); put("UUA", "L"); put("UUG", "L"); put("CUU", "L"); put("CUC", "L");
                put("CUA", "L"); put("CUG", "L"); put("AUU", "I"); put("AUC", "I"); put("AUA", "I"); put("AUG", "M");
                put("GUU", "V"); put("GUC", "V"); put("GUA", "V"); put("GUG", "V"); put("UCU", "S"); put("UCC", "S");
                put("UCA", "S"); put("UCG", "S"); put("CCU", "P"); put("CCC", "P"); put("CCA", "P"); put("CCG", "P");
                put("ACU", "T"); put("ACC", "T"); put("ACA", "T"); put("ACG", "T"); put("GCU", "A"); put("GCC", "A");
                put("GCA", "A"); put("GCG", "A"); put("UAU", "Y"); put("UAC", "Y"); put("UAA", "$"); put("UAG", "$");
                put("CAU", "H"); put("CAC", "H"); put("CAA", "Q"); put("CAG", "Q"); put("AAU", "N"); put("AAC", "N");
                put("AAA", "K"); put("AAG", "K"); put("GAU", "D"); put("GAC", "D"); put("GAA", "E"); put("GAG", "E");
                put("UGU", "C"); put("UGC", "C"); put("UGA", "$"); put("UGG", "W"); put("CGU", "R"); put("CGC", "R");
                put("CGA", "R"); put("CGG", "R"); put("AGU", "S"); put("AGC", "S"); put("AGA", "R"); put("AGG", "R");
                put("GGU", "G"); put("GGC", "G"); put("GGA", "G"); put("GGG", "G");
            }});
    private static final String START_CODON = "AUG";
    private static final String STOP_CODON = "$";
    private static final String FAILED_TO_FIND_COMPLEMENT_MRNA_NUCLEOTIDE = "Failed to find complement mrna nucleotide to ";
    private static final int TRIPLET_LENGTH = 3;

    private List<VariantContext> leftVcfList;
    private static int originalShift;
    private static Set<String> possibleFinalStrings = new HashSet<>();
    private static Set<String> proteins = new HashSet<>();
    private static PrintWriter resultsFile = null;

    public SequenceProcessor(List<VariantContext> leftVcfList, int originalShift, PrintWriter resultsFile) {
        this.leftVcfList = new ArrayList<>(leftVcfList);
        SequenceProcessor.originalShift = originalShift;
        SequenceProcessor.resultsFile = resultsFile;
    }

    public void process(String motherNucleotideString, String fatherNucleotideString, int motherStringShift,
                          int fatherStringShift, List<ExonPosition> motherExonPositions,
                        List<ExonPosition> fatherExonPositions)
            throws CloneNotSupportedException {
        StringBuilder motherNucleotideStringBuffer = new StringBuilder(motherNucleotideString);
        StringBuilder fatherNucleotideStringBuffer = new StringBuilder(fatherNucleotideString);

        while (leftVcfList.size() > 0) {
            VariantContext variantContext = leftVcfList.remove(0);

            if (variantContext.getGenotypes().get(0).isHet()) {
                if (variantContext.getGenotypes().get(0).isPhased()) {
                    String genotypeString = variantContext.getGenotypes().get(0).getGenotypeString(false);
                    if (genotypeString.indexOf("*") < genotypeString.indexOf("|")) {
                        fatherStringShift = applyMutation(variantContext, fatherNucleotideStringBuffer,
                                fatherStringShift, fatherExonPositions);
                    } else {
                        motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer,
                                motherStringShift, motherExonPositions);
                    }
                } else {
                    StringBuilder fatherNucleotideStringBufferForCopy = new StringBuilder(fatherNucleotideStringBuffer);
                    List<ExonPosition> fatherExonPositionsForCopy = fatherExonPositions.stream().map(ExonPosition::new).
                            collect(Collectors.toList());
                    int fatherStringShiftForCopy = applyMutation(variantContext, fatherNucleotideStringBufferForCopy,
                            fatherStringShift, fatherExonPositionsForCopy);

                    SequenceProcessor sequenceProcessorWithoutCurrentMutation = new SequenceProcessor(leftVcfList, originalShift, resultsFile);
                    sequenceProcessorWithoutCurrentMutation.process(motherNucleotideStringBuffer.toString(),
                            fatherNucleotideStringBufferForCopy.toString(), motherStringShift,
                            fatherStringShiftForCopy, motherExonPositions.stream().map(ExonPosition::new).
                                    collect(Collectors.toList()), fatherExonPositionsForCopy);

                    motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer, motherStringShift,
                            motherExonPositions);
                }
            } else {
                motherStringShift = applyMutation(variantContext, motherNucleotideStringBuffer, motherStringShift,
                        motherExonPositions);
                fatherStringShift = applyMutation(variantContext, fatherNucleotideStringBuffer, fatherStringShift,
                        fatherExonPositions);
            }
        }

        addNewResult(motherNucleotideStringBuffer, motherExonPositions);
        addNewResult(fatherNucleotideStringBuffer, fatherExonPositions);
    }

    public void addReferenceResult(String nucleotideString, List<ExonPosition> exonPositions) {
        addNewResult(new StringBuilder(nucleotideString), exonPositions);
    }

    private void addNewResult(StringBuilder nucleotideStringBuffer, List<ExonPosition> exonPositions) {
        if (!possibleFinalStrings.contains(nucleotideStringBuffer.toString())) {
            possibleFinalStrings.add(nucleotideStringBuffer.toString());

            Collections.sort(exonPositions, (ep1, ep2) -> ep1.getStartIndex() - ep2.getStartIndex());
            String nonInvertedMrnaString = generateNonInvertedMrnaString(nucleotideStringBuffer.toString());
            ExonPosition relativeStartCodonPosition = new ExonPosition();
            String invertedMrnaWithoutIntronsString = generateInvertedMrnaString(
                    new StringBuilder(generateMrnaWithoutIntronsString(nonInvertedMrnaString, exonPositions,
                            relativeStartCodonPosition)));
            String proteinString = generateProteinString(invertedMrnaWithoutIntronsString,
                    relativeStartCodonPosition.getStartIndex());

            resultsFile.println(nucleotideStringBuffer.toString());
            resultsFile.println(proteinString);
            resultsFile.println();

            if (!"".equals(proteinString)) {
                proteins.add(proteinString);
            }
        }
    }

    private int applyMutation(VariantContext variantContext, StringBuilder nucleotideStringBuilder, int shift,
                              List<ExonPosition> exonPositions) {
        nucleotideStringBuilder.replace(variantContext.getStart() - shift, variantContext.getEnd() - shift + 1,
                variantContext.getAlleles().get(1).getBaseString());

        int mutationShift = variantContext.getAlleles().get(0).getBaseString().length()
                - variantContext.getAlleles().get(1).getBaseString().length();

        for (ExonPosition exonPosition : exonPositions) {
            if (variantContext.getStart() < exonPosition.getStartIndex()) {
                exonPosition.setStartIndex(exonPosition.getStartIndex() - mutationShift);
            }
            if (variantContext.getStart() < exonPosition.getEndIndex()) {
                exonPosition.setEndIndex(exonPosition.getEndIndex() - mutationShift);
            }
        }

        return shift + mutationShift;
    }

    private String generateNonInvertedMrnaString(String nucleotideString) {
        StringBuilder mrnaStringBuilder = new StringBuilder();

        for (int i = 0; i < nucleotideString.length(); i++){
            char nucleotide = nucleotideString.charAt(i);
            Character complementMrnaNucleotide = COMPLEMENT_MRNA_NUCLEOTIDE_MAP.get(nucleotide);
            if (complementMrnaNucleotide != null) {
                mrnaStringBuilder.append(complementMrnaNucleotide);
            } else {
                LOGGER.error(FAILED_TO_FIND_COMPLEMENT_MRNA_NUCLEOTIDE + nucleotide);
            }
        }

        return mrnaStringBuilder.toString();
    }

    private String generateInvertedMrnaString(StringBuilder nonInvertedMrnaStringBuilder) {
        return nonInvertedMrnaStringBuilder.reverse().toString();
    }

    private String generateMrnaWithoutIntronsString(String nonInvertedMrnaString, List<ExonPosition> exonPositions,
                                                    ExonPosition relativeStartCodonPosition) {
        StringBuilder nucleotideStringBuilder = new StringBuilder();

        ExonPosition startCodon = new ExonPosition();
        ExonPosition exonWithStartCodon = new ExonPosition();
        Optional<ExonPosition> optionalStartCodon = exonPositions.stream().filter(exonPosition ->
                ExonType.START_CODON.equals(exonPosition.getExonType())).findFirst();
        if (optionalStartCodon.isPresent()) {
            startCodon.setExonPosition(optionalStartCodon.get());
            Optional<ExonPosition> optionalExonWithStartCodon = exonPositions.stream().filter(exonPosition ->
                    ExonType.EXON.equals(exonPosition.getExonType())
                            && exonPosition.getStartIndex() <= startCodon.getStartIndex()
                            && exonPosition.getEndIndex() >= startCodon.getEndIndex()).findFirst();
            if (optionalExonWithStartCodon.isPresent()) {
                exonWithStartCodon.setExonPosition(optionalExonWithStartCodon.get());
            }
        }

        exonPositions.stream().filter(exonPosition -> ExonType.EXON.equals(exonPosition.getExonType())).
                forEachOrdered(exonPosition -> {
                    if (exonWithStartCodon.equals(exonPosition)) {
                        relativeStartCodonPosition.setExonPosition(startCodon);
                        relativeStartCodonPosition.setStartIndex(nucleotideStringBuilder.length()
                                + startCodon.getStartIndex() - exonWithStartCodon.getStartIndex());
                        relativeStartCodonPosition.setEndIndex(nucleotideStringBuilder.length()
                                + startCodon.getEndIndex() - exonWithStartCodon.getStartIndex());
                    }

                    nucleotideStringBuilder.append(nonInvertedMrnaString.substring(exonPosition.getStartIndex()
                                    - originalShift, exonPosition.getEndIndex() - originalShift + 1));
                });

        return generateInvertedMrnaString(nucleotideStringBuilder);
    }

    private String generateProteinString(String invertedMrnaWithoutIntrons, int startCodonIndex) {
        StringBuilder proteinStringBuilder = new StringBuilder();
        StringBuilder mrnaFromStartCodonStringBuilder = new StringBuilder(invertedMrnaWithoutIntrons.substring(
                0, startCodonIndex + TRIPLET_LENGTH));
        mrnaFromStartCodonStringBuilder = mrnaFromStartCodonStringBuilder.reverse();

        if (!START_CODON.equals(mrnaFromStartCodonStringBuilder.substring(0, TRIPLET_LENGTH))) {
            return "";
        }

        for (int i = 0; i + 2 < mrnaFromStartCodonStringBuilder.length(); i += TRIPLET_LENGTH) {
            String aminoAcid = NUCLEOTIDE_TRIPLET_AMINO_ACID_MAP.get(mrnaFromStartCodonStringBuilder.
                    substring(i, i + TRIPLET_LENGTH));
            if (STOP_CODON.equals(aminoAcid)) {
                return proteinStringBuilder.toString();
            }
            proteinStringBuilder.append(aminoAcid);
        }

        return proteinStringBuilder.toString();
    }

    public static Set<String> getProteins() {
        return proteins;
    }
}
