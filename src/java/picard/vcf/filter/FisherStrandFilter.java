package picard.vcf.filter;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Filters records based on the phred scaled p-value from the Fisher Strand test stored in
 * the FS attribute.
 *
 * @author tfennell
 */
public class FisherStrandFilter implements VariantFilter {
    private double maxPhredScalePValue;

    public FisherStrandFilter(final double maxPhredScalePValue) {
        this.maxPhredScalePValue= maxPhredScalePValue;
    }

    @Override
    public String filter(final VariantContext ctx) {
        final double fs = ctx.getAttributeAsDouble("FS", 0);
        return (fs > maxPhredScalePValue) ? "StrandBias" : null;
    }
}
