package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/** Focused native Mutect2 coverage of the shared phasing path on fabricated reads. */
public class Mutect2SyntheticPhaseIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testMultiallelicDeletionAndSupportedSnpPhase() {
        final File bam = new File(toolsTestDir, "walkers/GenotypeGVCFs/native-extra-alt.bam");
        final File reference = new File(toolsTestDir, "walkers/GenotypeGVCFs/native-extra-alt.fa");
        final File output = createTempFile("syntheticSomaticDeletionPhase", ".vcf");
        runCommandLine(new ArgumentsBuilder().addInput(bam).addReference(reference)
                .addInterval(new SimpleInterval("fabricated:40-120"))
                .add("tumor-sample", "synthetic")
                .add(M2ArgumentCollection.EMISSION_LOG_SHORT_NAME, 0).addOutput(output), "Mutect2");

        final Map<Integer, VariantContext> calls = VariantContextTestUtils.getVariantContexts(output).stream()
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        Assert.assertEquals(calls.keySet(), Set.of(55, 61, 65));
        final VariantContext upstream = calls.get(55);
        final VariantContext parent = calls.get(61);
        final VariantContext downstream = calls.get(65);
        Assert.assertEquals(parent.getReference().getBaseString(), "GAAGTCC");
        Assert.assertEquals(parent.getAlternateAlleles(), List.of(Allele.create("G"), Allele.create("GC")));
        Assert.assertEquals(parent.getGenotype(0).getAlleles(), parent.getAlleles());
        Assert.assertFalse(parent.getGenotype(0).isPhased());
        Assert.assertFalse(parent.getGenotype(0).hasExtendedAttribute(VCFConstants.PHASE_SET_KEY));
        Assert.assertEquals(downstream.getAlternateAlleles(), List.of(Allele.create("G")));
        Assert.assertFalse(downstream.getAlternateAlleles().contains(Allele.SPAN_DEL));
        Assert.assertEquals(downstream.getGenotype(0).getAlleles(),
                List.of(downstream.getAlternateAllele(0), downstream.getReference()));
        Assert.assertEquals(upstream.getGenotype(0).getAlleles(),
                List.of(upstream.getReference(), upstream.getAlternateAllele(0)));
        Assert.assertTrue(upstream.getGenotype(0).isPhased());
        Assert.assertTrue(downstream.getGenotype(0).isPhased());
        Assert.assertEquals(String.valueOf(upstream.getGenotype(0).getExtendedAttribute(VCFConstants.PHASE_SET_KEY)), "55");
        Assert.assertEquals(upstream.getGenotype(0).getExtendedAttribute(VCFConstants.PHASE_SET_KEY),
                downstream.getGenotype(0).getExtendedAttribute(VCFConstants.PHASE_SET_KEY));
    }
}
