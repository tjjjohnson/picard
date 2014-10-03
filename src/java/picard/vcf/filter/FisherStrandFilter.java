package picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import java.util.List;

/**
 * Filters records based on the phred scaled p-value from the Fisher Strand test stored in
 * the FS attribute.
 *
 * @author tfennell
 */
public class FisherStrandFilter implements VariantFilter {
    private final double maxPhredScalePValue;

    public FisherStrandFilter(final double maxPhredScalePValue) {
        this.maxPhredScalePValue= maxPhredScalePValue;
    }

    @Override
    public List<VCFFilterHeaderLine> headerLines() {
        return CollectionUtil.makeList(new VCFFilterHeaderLine("StrandBias", "Site exhibits excessive allele/strand correlation."));
    }

    @Override
    public String filter(final VariantContext ctx) {
        final double fs = ctx.getAttributeAsDouble("FS", 0);
        return (fs > maxPhredScalePValue) ? "StrandBias" : null;
    }
}
