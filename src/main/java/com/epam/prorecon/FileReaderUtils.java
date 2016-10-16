package com.epam.prorecon;

import com.epam.prorecon.entity.ExonPosition;
import com.epam.prorecon.entity.ExonType;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class FileReaderUtils {
    private static final Logger LOGGER = LoggerFactory.getLogger(FileReaderUtils.class);
    private static final String UTF_8 = "UTF-8";

    private FileReaderUtils() {
    }

    public static String readSequenceFromFastaFile(String fileUrl, String chromosomeName, Integer beginIndex,
                                                   Integer endIndex) throws FileNotFoundException {
        final File refFile = new File(fileUrl);
        IndexedFastaSequenceFile sequenceFile = new IndexedFastaSequenceFile(refFile);

        String nucleotideString = sequenceFile.getSubsequenceAt(chromosomeName, beginIndex, endIndex).getBaseString();
        CloserUtil.close(sequenceFile);
        return nucleotideString;
    }

    public static List<VariantContext> readVariantContextsFromVcfFile(String fileUrl, String chromosomeName,
                                                                      Integer beginIndex, Integer endIndex) {
        final File vcfFile = new File(fileUrl);
        VCFFileReader vcfFileReader = new VCFFileReader(vcfFile, false);

        List<VariantContext> variantContexts = vcfFileReader.query(chromosomeName, beginIndex, endIndex).toList();
        CloserUtil.close(vcfFile);
        return variantContexts;
    }

    public static List<ExonPosition> readGffFile(String fileUrl) {
        List<ExonPosition> exonPositions = new ArrayList<>();
        String line;
        try (
                InputStream fis = new FileInputStream(fileUrl);
                InputStreamReader isr = new InputStreamReader(fis, Charset.forName(UTF_8));
                BufferedReader br = new BufferedReader(isr)
        ) {
            while ((line = br.readLine()) != null) {
                Matcher m = Pattern.compile("(\\w+)\\s+\\w+\\s+(\\w+)\\s+(\\d+)\\s+(\\d+)").matcher(line);

                while (m.find()) {
                    ExonType exonType = ExonType.asExonType(m.group(2).replace("5", "five_").replace("3", "three_"));
                    if (exonType != null) {
                        exonPositions.add(new ExonPosition(m.group(1), exonType, Integer.parseInt(m.group(3)),
                                Integer.parseInt(m.group(4))));
                        if (exonType.equals(ExonType.THREE_UTR)) {
                            return exonPositions;
                        }
                    }
                }
            }
        } catch (IOException e) {
            LOGGER.error(e.getMessage());
        }

        return exonPositions;
    }
}
