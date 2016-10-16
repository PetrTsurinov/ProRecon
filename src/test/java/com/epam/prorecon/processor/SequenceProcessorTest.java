package com.epam.prorecon.processor;

import com.epam.prorecon.FileReaderUtils;
import com.epam.prorecon.entity.ExonPosition;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.Assert;
import org.junit.Test;

import java.io.*;
import java.net.URL;
import java.util.List;
import java.util.stream.Collectors;

public class SequenceProcessorTest {
    @Test
    public void nucliotideStringAfterApplyingVcfShouldBeEqualToPrecalculatedConstant() throws IOException, CloneNotSupportedException {
        URL fastaFileUrl = this.getClass().getClassLoader().getResource("dmel-all-chromosome-r606.fasta");
        Assert.assertNotNull(fastaFileUrl);
        String fastaFileSubSequence =  FileReaderUtils.readSequenceFromFastaFile(fastaFileUrl.getPath(), "X", 12584385, 12592193);

        URL vcfFileUrl = this.getClass().getClassLoader().getResource("agnX1.model.2.snp-indels.vcf");
        Assert.assertNotNull(vcfFileUrl);
        File vcfFile = new File(vcfFileUrl.getPath());

        URL gtfFileUrl = this.getClass().getClassLoader().getResource("dmel-all-r6.06.LIMK1.gtf");
        Assert.assertNotNull(gtfFileUrl);

        File vcfIndexFile = new File(vcfFileUrl.getPath() + ".Idx");
        Index idx = IndexFactory.createIndex(vcfFile, new VCFCodec(), IndexFactory.IndexType.LINEAR);

        LittleEndianOutputStream stream;
        stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(vcfIndexFile)));
        idx.write(stream);
        stream.close();

        CloserUtil.close(vcfFile);
        CloserUtil.close(vcfIndexFile);

        PrintWriter totalFile = new PrintWriter("total.txt");
        List<ExonPosition> exonPositions = FileReaderUtils.readGffFile(gtfFileUrl.getPath());
        SequenceProcessor sequenceProcessor = new SequenceProcessor(FileReaderUtils.readVariantContextsFromVcfFile(
                vcfFileUrl.getPath(), "X", 12584385, 12592193), 12584385, totalFile);

        sequenceProcessor.addReferenceResult(fastaFileSubSequence, exonPositions);
        sequenceProcessor.process(fastaFileSubSequence, fastaFileSubSequence, 12584385, 12584385, exonPositions,
                exonPositions.stream().map(ExonPosition::new).collect(Collectors.toList()));
        totalFile.close();

        Assert.assertEquals(19, SequenceProcessor.getProteins().size());
    }
}
