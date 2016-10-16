package com.epam.prorecon;

import com.epam.prorecon.entity.ExonPosition;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.junit.Assert;
import org.junit.Test;

import java.io.*;
import java.net.URL;
import java.util.List;

public class FileReaderUtilsTest {
    private final static String FIRST_FASTA_SEQUANCE_BEGIN =
            "TTAAAATTCCACAATAAAGACACTAGAAAGTTTAAAAGTTTTTTTTTAATTTCCTCTTTTATTAACAAAAAA" +
                    "TGTTTTAAGATTAAAATCAAAAGAGAATCTACGAACTACGTTTATACACAAGC";

    @Test
    public void readSequenceFromFastaFileTest() throws FileNotFoundException {
        URL fastaFileUrl = this.getClass().getClassLoader().getResource("dmel-all-chromosome-r606.fasta");
        Assert.assertNotNull(fastaFileUrl);
        String fastaFileSubSequence = FileReaderUtils.readSequenceFromFastaFile(
                fastaFileUrl.getPath(), "X", 12584276, 12584276 + FIRST_FASTA_SEQUANCE_BEGIN.length() - 1);

        Assert.assertEquals(FIRST_FASTA_SEQUANCE_BEGIN, fastaFileSubSequence);
    }

    @Test
    public void readVariantContextsFromVcfFile() throws IOException {
        URL vcfFileUrl = this.getClass().getClassLoader().getResource("agnX1.model.2.snp-indels.vcf");
        Assert.assertNotNull(vcfFileUrl);
        File vcfFile = new File(vcfFileUrl.getPath());

        File vcfIndexFile = new File(vcfFileUrl.getPath() + ".Idx");
        Index idx = IndexFactory.createIndex(vcfFile, new VCFCodec(), IndexFactory.IndexType.LINEAR);

        LittleEndianOutputStream stream;
        stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(vcfIndexFile)));
        idx.write(stream);
        stream.close();

        CloserUtil.close(vcfFile);
        CloserUtil.close(vcfIndexFile);

        List<VariantContext> variantContexts = FileReaderUtils.readVariantContextsFromVcfFile(
                vcfFileUrl.getPath(), "X", 12585170, 12585230);

        Assert.assertNotNull(variantContexts);
        Assert.assertEquals(2, variantContexts.size());
        Assert.assertNotNull(variantContexts.get(0));
        List<Allele> firstVcAlleles = variantContexts.get(0).getAlleles();
        Assert.assertEquals(2, firstVcAlleles.size());
        Assert.assertEquals("CA", firstVcAlleles.get(0).getBaseString());
        Assert.assertEquals("C", firstVcAlleles.get(1).getBaseString());
        List<Allele> secondVcAlleles = variantContexts.get(1).getAlleles();
        Assert.assertEquals(2, secondVcAlleles.size());
        Assert.assertEquals("CT", secondVcAlleles.get(0).getBaseString());
        Assert.assertEquals("C", secondVcAlleles.get(1).getBaseString());
    }

    @Test
    public void readGtfFileTest() throws FileNotFoundException {
        URL gtfFileUrl = this.getClass().getClassLoader().getResource("dmel-all-r6.06.LIMK1.gtf");
        Assert.assertNotNull(gtfFileUrl);
        List<ExonPosition> exonPositions = FileReaderUtils.readGffFile(gtfFileUrl.getPath());

        Assert.assertNotNull(exonPositions);
    }
}
